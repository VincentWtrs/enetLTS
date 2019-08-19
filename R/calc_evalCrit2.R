# calc_evalCrit() Allowing nfold = 1 FOR INFORMATION CRITERION-BASED APPROACHES
calc_evalCrit2 <- function(rowind, combis_ind, alphas, lambdas, 
                           index, xx, yy, nfold, repl, family, ic_type = NULL) {
  # family argument defined as well, because it was defined within cv.enetLTS and calc_evalCrit just used it there as well
  
  ## NEW: Handling information criterion case ## NEW(2) moved the nfold and repl correcting part to top-layer enetLTS() function
  if(!is.null(ic_type)){
    ic <- TRUE 
  } else if (is.null(ic_type)) {
    ic <- FALSE
  }
  
  # EXtracting the lambda and alpha currently under analysis
  i <- combis_ind[rowind, 1]
  j <- combis_ind[rowind, 2]
  lambda <- lambdas[i]
  alpha <- alphas[j]
  
  # Printing which alpha-lambda combination is currently used
  if(isFALSE(ic)){
    print(paste("Tuning with cross-validation approach, currently fitting for: ", alpha, " and lambda :", lambda), sep = "")
  } else if (isTRUE(ic)) {
    print(paste("Tuning with information criterion approach, currently fitting for: ", alpha, " and lambda :", lambda), sep = "")
  }
  
  ## Index handling
  # No index: use all data
  if (is.null(index)) {
    x <- xx
    y <- yy
    # Yes index: use the subsetted data by the given indices
  } else {
    x <- xx[index[, i, j], ]
    y <- yy[index[, i, j]]
  }
  
  ### LOOPING OVER REPS AND FOLDS
  ## For each replication
  evalCritl <- rep(NA, repl) # Initiating
  for (l in 1:repl) {
    
    ## NEW: if(ic == FALSE): Keep old functionality
    if (isFALSE(ic)) {
      # Binomial Case
      if (family == "binomial") {
        folds0 <- cvTools:::cvFolds(length(y[y == 0]), K = nfold, R = 1, type = "random")
        folds1 <- cvTools:::cvFolds(length(y[y == 1]), K = nfold, R = 1, type = "random")
        loss0 <- rep(NA, sum(y == 0))
        loss1 <- rep(NA, sum(y == 1))
        
        # Gaussian Case
      } else if (family == "gaussian") {
        folds <- cvTools:::cvFolds(length(y), K = nfold, R = 1, type = "random")
        loss <- rep(NA, nrow(x))
      }
      
      #### Cross Validation (CV) Loop
      for (f in 1:nfold) {
        ### Creating folds
        ## Binomial Case
        if (family == "binomial") {
          
          # 0-outcomes
          xtrain0 <- x[y == 0, ][folds0$subsets[folds0$which != f, 1], ]
          ytrain0 <- y[y == 0][folds0$subsets[folds0$which != f, 1]]
          xtest0 <- x[y == 0, ][folds0$subsets[folds0$which == f, 1], ]
          ytest0 <- y[y == 0][folds0$subsets[folds0$which == f, 1]]
          
          # 1-outcomes
          xtrain1 <- x[y == 1, ][folds1$subsets[folds1$which != f, 1], ]
          ytrain1 <- y[y == 1][folds1$subsets[folds1$which !=  f, 1]]
          xtest1 <- x[y == 1, ][folds1$subsets[folds1$which == f, 1], ]
          ytest1 <- y[y == 1][folds1$subsets[folds1$which == f, 1]]
          
          # Combining
          xtrain <- rbind(xtrain0, xtrain1)
          ytrain <- c(ytrain0, ytrain1)
          xtest <- rbind(xtest0, xtest1)
          ytest <- c(ytest0, ytest1)
          
          ## Gaussian Case  
        }  else if (family == "gaussian") {
          xtrain <- x[folds$subsets[folds$which != f, 1], ]
          ytrain <- y[folds$subsets[folds$which != f, 1]]
          xtest <- x[folds$subsets[folds$which == f, 1], ]
          ytest <- y[folds$subsets[folds$which == f, 1]]
        }
        
        ### Fitting elastic net for each fold
        res <- tryCatch({
          hpen <- length(ytrain)
          if(family == "binomial"){
            trainmod <- glmnet(x = xtrain, 
                               y = ytrain, 
                               family = family, 
                               lambda = lambda/hpen, # hpen !=h because of folds!
                               alpha = alpha, 
                               standardize = FALSE, # Because already robustly standardized!
                               intercept = TRUE,  # WITH INTERCEPT BECAUSE BINOMIAL
                               maxit = 10^6) # Allowing for additional extra iterations (default: 10^5) in case of nonconvergence
          } else if(family == "gaussian"){
            trainmod <- glmnet(x = xtrain, 
                               y = ytrain, 
                               family = family, 
                               lambda = lambda/hpen, # hpen !=h because of folds!
                               alpha = alpha, 
                               standardize = FALSE, # Because already robustly standardized!
                               intercept = FALSE) # WITHOUT BECAUSE GAUSSIAN, WILL GO THROUGH ORIGIN OR CLOSE
          }
          
        }, error = function(err) {
          error <- TRUE
          return(error)
        })
        if (is.logical(res)) {
          print(paste("CV broke off for alpha=", alpha, "and lambda=", lambda))
        }
        else {
          trainmod <- res
          if (family == "binomial") {
            loss0[folds0$which == f] <- -(ytest0 * xtest0 %*% matrix(trainmod$beta)) + log(1 + exp(xtest0 %*% matrix(trainmod$beta))) # ORIGINAL
            loss1[folds1$which == f] <- -(ytest1 * xtest1 %*% matrix(trainmod$beta)) + log(1 + exp(xtest1 %*% matrix(trainmod$beta))) # ORIGINAL
          }
          else if (family == "gaussian") 
            loss[folds$which == f] <- ytest - xtest %*% 
              matrix(trainmod$beta)
        }
      }
      if (family == "binomial") {
        loss <- c(loss0, loss1)
        evalCritl[l] <- mean(loss, na.rm = TRUE)
      }
      else if (family == "gaussian") {
        evalCritl[l] <- sqrt(mean(loss^2))
      }
    } # END OF if(ic == FALSE)
    
    ########## NEW: if ic is TRUE ########## 
    if (isTRUE(ic)) {
      loss <- rep(NA, nrow(x)) # nrow(x) = h
      xtrain <- x # Note: there is only a training set here
      ytrain <- y # Note: there is only a training set here
      
      # Fitting (within an error catching structure)
      res <- tryCatch({
        hpen <- length(ytrain) # Sample size WITHIN THE FOLD (< h)
        
        if (family == "binomial"){
          trainmod <- glmnet(xtrain, 
                             ytrain, 
                             family, 
                             alpha = alpha, 
                             lambda = lambda/hpen, # Note: hpen != h
                             standardize = FALSE, # NOT NEEDED BECAUSE ALREADY STANDARDIZED
                             intercept = TRUE, # BECAUSE BINOMIAL
                             maxit = 10^6) 
          
        } else if (family == "gaussian") {
          trainmod <- glmnet(xtrain, 
                             ytrain, 
                             family, 
                             alpha = alpha, 
                             lambda = lambda/hpen, # Note: hpen != h
                             standardize = FALSE, # NOT NEEDED BECAUSE ALREADY STANDARDIZED
                             intercept = FALSE) # BECAUSE STANDARDIZATION WILL CAUSE IT TO GO THROUGH ORIGIN OR CLOSE
        }
        
      }, error = function(err) {
        error <- TRUE
        return(error)
      })
      if (is.logical(res)) {
        print(paste("Fitting broke off for alpha =", alpha, "and lambda =", lambda))
      }
      else {
        trainmod <- res
        # Binomial Case
        if (family == "binomial") {
          # NEW
          loss <- -(y * xtrain %*% matrix(trainmod$beta)) + log(1 + exp(xtrain %*% matrix(trainmod$beta))) # Loss (Neg. loglik, length = h)
        }
        
        # Gaussian Case
        else if (family == "gaussian") 
          loss <- ytrain - xtrain %*% matrix(trainmod$beta)
      }
      
      ## Calculating evaluation criteria
      # Binomial Case
      if (family == "binomial") {
        evalCritl[l] <- mean(loss, na.rm = TRUE)
        if(ic == TRUE){
          evalCritl[l] <- 2 * mean(loss, na.rm = TRUE) + ic_penalty(model = trainmod, 
                                                                    type = ic_type, 
                                                                    X = xtrain, 
                                                                    alpha = alpha,
                                                                    intercept = TRUE) # BECAUSE BINOMIAL ALWAYS WITH INTERCEPT
          # Note: the intercept setting needs to be the same as used to fit the respective model that is given to ic_penalty
        }
      } # Gaussian Case
      else if (family == "gaussian") {
        evalCritl[l] <- sqrt(mean(loss^2))
        if(ic == TRUE) {
          if (length(ic_type) == 1) {
            evalCritl[l] <- 2 * evalCritl[l] + ic_penalty(model = trainmod, 
                                                          type = ic_type, 
                                                          X = xtrain, 
                                                          alpha = alpha,
                                                          intercept = FALSE) # Minus loss + Penalty
            print(evalCritl)
          } else if (length(ic_type) > 1) {
            evalCritl <- rep(NA, length = length(ic_type))
            for (i in 1:length(ic_type)) {
              evalCritl[i] <-  evalCritl[l] + ic_penalty(model = trainmod, 
                                                         type = ic_type, 
                                                         X = xtrain, 
                                                         alpha = alpha,
                                                         intercept = FALSE) # Minus loss + Penalty
            }
          }
        }
      }
    } 
    # END OF if(ic == TRUE)
  } # END OF REPL. LOOP
  return(list(lambda_ind = i, alpha_ind = j, evalCritl = evalCritl))
}