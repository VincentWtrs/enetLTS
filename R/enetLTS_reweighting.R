enetLTS_raw_fit <- function(xx, yy, family, indexbest, alphabest, lambdabest, intercept){
  
  # In this part we will get a final fit. Common practice in machine learning, after the optimal hyperparameters have been found
  # ... to take this hyperparameter and fit it to all data. Here, this case, it is slightly more complex: there is (i) the tuned model
  # ... that is "non-reweighted", for which a final fit can be obtained. Secondly (ii) there is the reweighted tuned model for which
  # ... a final fit can be obtained as well. The issue is that not all data is used but only the 'supposed' outlier-free set, of course.
  # ... The reweighting step also complicates things even further since the amount of outliers will be reduced.
  
  
  #### FINAL FITTING OF THE MODEL USING THE OPTIMAL PARAMETERS (NON-REWEIGHTED FIT)
  
  ### STEPS FOR SCAL: TRUE
  ## STEP 1. Scaling (if scal TRUE) using the best indices found (related to best hyperparameter), hence nonrobust (because outlier-free world)
  ## STEP 2. Fitting of ordinary elastic net (glmnet) on all non-outlying data with best hyperparams and scaled outlier-free dataset
  ## STEP 3. Handling of the results from the final (non-reweighted) fit.
  # STEP 3a. Intercept handling of (non-reweighted) fit
  # STEP 3b. Coefficients handling of (non-reweighred) fit
  # STEP 3c. Residuals handling of (non-reweighted) fit
  #---------------- REWEIGHTING ----------------------#
  ## STEP 4: Calculating raw weights (raw.wt)
  ## STEP 5: Scaling (if scal TRUE) using the best REWEIGHTED indices found found (related to best hyperparamters), hence nonrobust but potentially less conservative
  
  ## STEP 7
  ## STEP 8b. INTERCEPT HANDLING OF REWEIGHTED FIT
  ## STEP 8b. COEFFICIENT HANDLING OF REWEIGHTED FIT
  
  
  
  ### STEP 1. Scaling of data
  if (isTRUE(scal)) {
    scl <- prepara(x = xx, 
                   y = yy, 
                   family = family, 
                   index = indexbest, # Indices provided s.t. the normalization should happen on the outlier-free set (=='outlier-free world')
                   robu = 0) # Nonrobust because now in outlier-free world, as provided by the indices above
    xs <- scl$xnor # Normalized X
    ys <- scl$ycen # Centered y (normally, but since family is passed on, no centering is done for Binomial!)
    
    ### STEP 2: Final fitting elastic net with optimal hyperparameters and outlier-free set (non-reweighted)
    ## STEP 2: Final fitting: Case Binomial
    if (family == "binomial") {
      fit <- glmnet(x = xs[indexbest, ],
                    y = ys[indexbest, ],
                    family = "binomial",
                    alpha = alphabest,  # With the tuned alpha
                    lambda = lambdabest,  # ... and the tuned lambda
                    standardize = FALSE, # Because already done in the prepara case (we are providing y standardized and x standardized)
                    intercept = TRUE) # Because only for linear models the standardization makes the regression line go through the origin, not true for GLMs in general
      
    ## STEP 2: Final fitting: Case Gaussian
    } else if (family == "gaussian") {
      fit <- glmnet(x = xs[indexbest, ], 
                    y = ys[indexbest, ], 
                    family = "gaussian", 
                    alpha = alphabest, 
                    lambda = lambdabest, 
                    standardize = FALSE, # Because already done
                    intercept = FALSE)  # FALSE becasue for standardized data, a linear model goes through the origin!
    }
    
    ### STEP 3: HANDLING RESULTS OF THE FINAL (NON-REWEIGHTED) FIT
    ## STEP 3a. Intercept Handling
    # STEP 3a. Intercept handling: Case Binomial
    if (family == "binomial") {
      if (isFALSE(intercept)) {
        a00 <- 0 # NEW: Changed this from weird a00 <- if (intercept == FALSE) {0} style of notation, this is clearer, does the same
      } else {
        a00 <- fit$a0 # NEW: I think this is the way (a0 is the way a glmnet object can have its intercept accessed)
        #a00 <- drop(fit$a0 - as.vector(as.matrix(fit$beta)) %*% (scl$mux / scl$sigx)) # NEW: REMOVED SINCE THIS WOULD ONLY! MAKE SENSE IN LINEAR CASES
        # NOTE: This only works if the model is fitted with a single lambda otherwise $a0 will give a vector of intercepts
      }
      
      # STEP 3b. Coefficients handling of (non-reweighted) fit
      raw.coefficients <- drop(as.matrix(fit$beta) / scl$sigx) # Getting original coefs (w.r.t. unstandardized data)
      # NOTE: Formula above holds for all GLMs due to invariance principle. Note that this must be done since the data provided was standardized before the function! I think here we could just run glmnet as-is as well...
      beta_with_int <- drop(as.matrix((coef(fit, s = lambdabest)))) # Adding s = lambdabest as redundancy, if no lambdabest would be supplied to glmnet, it would be able to display it for all lambdas
      # Note calling coef on a glmnet object will give those dcg sparse matrices. By forcing to matrix and dropping unnecessary dimensions we get a named vector
      
      # STEP 3c. Residuals handling of (non-reweighted) fit
      raw.residuals <- -(ys * cbind(1, xs) %*% beta_with_int) + log(1 + exp(cbind(1, xs) %*% beta_with_int)) # NEW: cbind(1, xs) AND beta_with_int TO TAKE INTERCEPT INTO ACCOUNT
      #raw.residuals <- -(ys * xs %*% as.matrix(fit$beta)) + log(1 + exp(xs %*% as.matrix(fit$beta))) # OLD FORMULA, NOT CORRECT FOR BINOMIAL SINCE IT DOESN'T TAKE INTERCEPT INTO ACCOUNT
      
      ##################### REWEIGHTING STEP #####################
      
      ## STEP 4: CALCULATING RAW WEIGHTS
      raw.wt <- weight.binomial(x = xx, 
                                y = yy, 
                                beta = c(a00, raw.coefficients),  # Or beta_with_int
                                intercept = intercept, 
                                del = del) # del = 0.125 by standard
      
      ## STEP 5: Scaling data, now using reweighted weights (based on STEP 4.)
      sclw <- prepara(x = xx, 
                      y = yy, 
                      family = family, 
                      index = which(raw.wt == 1), # Only for those which are outlier-free
                      robu = 0)  # Indeed non-robust since this should be outlier-free as well
      xss <- sclw$xnor # Extracting standardized X
      yss <- sclw$ycen # Extracting centered y (but for Binomial no centering is actually performed)
      # NOTE: This gives us a new dataset, for which we -in essence- can tune again
      
      ## STEP 6a: FITTING ALL HYPERPARAMETERS THEIR MODELS ()
      # STEP 6a: Fitting all hyperparameters their models if no lambdaw is given (full regularization path for alphas)
      if (missing(lambdaw)) {
        if (family == "binomial"){
          reweighted_cv <- cva.glmnet(x = xss[which(raw.wt == 1), ],
                                      y = yss[which(raw.wt == 1)],
                                      family = "binomial",
                                      nfolds = 10, # TODO: CHOOSE
                                      type.measure = "deviance", # This should be passable to lower lying functions such as cv.glmnet
                                      standardize = FALSE, # Since data already standardized!
                                      intercept = TRUE) # True for binomial!
          
          ## OLD CODE NOT TUNING FOR ALPHA
          #reweighted_cv <- cv.glmnet(x = xss[which(raw.wt == 1), ], # NEW: changed name to lambdaw -> lambdaw_fit
          #                         y = yss[which(raw.wt == 1)], 
          #                         family = family, 
          #                         nfolds = 5, 
          #                         alpha = alphabest, 
          #                         standardize = FALSE, 
          #                         intercept = TRUE, # NEW: changed this to true
          #                         type.measure = "deviance")
        } else if (family == "gaussian") {
          reweighted_cv <- cva.glmnet(x = xss[which(raw.wt == 1), ],
                                      y = yss[which(raw.wt == 1)],
                                      family = "gaussian",
                                      nfolds = 10,
                                      typ.measure = "mse", # This should be passable to lower lying functions such as cv.glmnet. Quid: is deviance basically the same in this case?
                                      standardize = FALSE,
                                      intercept = FALSE) # Since for Gaussian, this actually holds when having standardized the data
        }
        
        # TODO Maybe extract the used lambdaw also
        
      # STEP 6a: 'Tuning' hyperparameters with if single lambdaw is given (UNLIKELY, i.e. the best lambdaw is the only one given!)
      } else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw_best <- lambdaw
        
      # STEP 6a. Tuning hyperparamters if vector of lambdaw is given (UNLIKELY)
      } else if (!missing(lambdaw) & length(lambdaw) > 1) { # Multiple lambdaw given
        if (family == "binomial"){
          reweighted_cv <- cva.glmnet(x = xss[which(raw.wt == 1), ],  # With the new reweighted dataset
                                      y = yss[which(raw.wt == 1)],  # With the new reweighted dataset
                                      family = "binomial", 
                                      lambda = lambdaw, # Passing on lambdaw vector given by user
                                      nfolds = 10, # TODO CHOOSE
                                      type.measure = "deviance",  # For binomial models negative loglikelihood
                                      standardize = FALSE, 
                                      intercept = TRUE)
        } else if (family == "gaussian") {
          reweighted_cv <- cva.glmnet(x = xss[which(raw.wt == 1), ],
                                      y = yss[which(raw.wt == 1)],
                                      family = "gaussian",
                                      lambda = lambdaw, # Passing on lambdaw vector given by user
                                      nfolds = 10, # TODO CHOOSE
                                      type.measure = "mse",  # Usually for Gaussian models
                                      standardize = FALSE,
                                      intercept = FALSE)
        }
      }
      
      ## STEP 6b. CHOOSING OPTIMAL ALPHA (JUST MIN)
      alphaw_index <- which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm))) # Extracting index (bit weird but OK)
      alphaw_best <- reweighted_cv$alpha[alphaw_index] # Actual value of the best alpha!
      
      ## STEP 6c. CHOOSING OPTIMAL LAMBDA WITHIN OPTIMAL ALPHA # TODO: MIGHT DO THIS MORE INDEPENDENTLY
      if (type_lambdaw == "min") {
        lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.min # NEW!
        #lambdaw <- lambdaw_fit$lambda.min
        #lambdaw <- reweighted_cv$lambda[which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm)))]  # TODO (TEMP) Correct, this is to see if it works
      } else if (type_lambdaw == "1se") {
        lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.1se
        #lambdaw <- lambdaw_fit$lambda.1se
        #lambdaw <- reweighted_cv$lambda[which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm)))]  # TODO (TEMP) Correct, this is to see if it works
      }
      
      ## STEP 7: FITTING FINAL MODEL USING OPTIMAL HYPERPARAMETERS ON REWEIGHTED DATA
      print("alphabest:") # TODO REMOVE
      print(alphabest) # TODO REMVOE
      print("lambdaw:") # TODO REMOVE
      print(lambdaw) # TODO REMOVE
      
      # STEP 7: Binomial fitting final model using optimal hyperparameters on reweighted data
      if (family == "binomial"){
        fitw <- glmnet(x = xss[which(raw.wt == 1), ], 
                       y = yss[which(raw.wt == 1)], 
                       family = family, 
                       alpha = alphabest, 
                       lambda = lambdaw, 
                       standardize = FALSE, 
                       intercept = TRUE) # NEW: changed this to TRUE for binomial because standardization doesn't work that way for binomial
      }
      # STEP 7: Gaussian fitting final model using optimal hyperparameter on reweighted data
      if (family == "gaussian"){
        fitw <- glmnet(x = xss[which(raw.wt == 1), ], 
                       y = yss[which(raw.wt == 1)], 
                       family = family, 
                       alpha = alphabest, 
                       lambda = lambdaw, 
                       standardize = FALSE, 
                       intercept = FALSE) # Keeping at FALSE for Gaussian since using standardize data
      }
      
      # TODO REMOVE
      print("I am past fitting the final reweighted fit object")
      
      ## STEP 8a. INTERCEPT HANDLING OF REWEIGHTED FIT
      if (isFALSE(intercept)) {
        a0 <- 0
      } else { # intercept TRUE
        if (family == "binomial") {
          a0 <- fitw$a0  # Just get it from the model object
        } else if (family == "gaussian") {
          a0 <- drop(fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% (sclw$mux/sclw$sigx)) # For Gaussian it holds
        }
      }
      
      ## STEP 8b. COEFFICIENT HANDLING OF REWEIGHTED FIT
      coefficients <- drop(as.matrix(fitw$beta)/sclw$sigx)
      
      # TODO REMOVE SOME PRINT STATEMENTS
      print("Structure of coefficients:") # TODO REMOVE
      print(str(coefficients)) # # TODO REMOVE
      print("Structure of a0:") # # TODO REMOVE
      print(str(a0)) # TODO REMOVE
      print("Calculated coefs")  # TODO REMOVE
      
      ## STEP 8c. CALCULATING WEIGHTS
      wgt <- weight.binomial(x = xx, 
                             y = yy, 
                             beta = c(a0, coefficients), 
                             intercept = intercept, 
                             del = del)
      
      reweighted.residuals <- -(yy * cbind(1, xx) %*% c(a0, coefficients)) + log(1 + exp(cbind(1, xx) %*% c(a0, coefficients)))
      print("Calculated reweighed resids")  # TODO REMoVE
      
      # Case: Gaussian
    } else if (family == "gaussian") {
      a00 <- if (intercept == FALSE) 
        0
      else drop(scl$muy + fit$a0 - as.vector(as.matrix(fit$beta)) %*% (scl$mux/scl$sigx))
      
      raw.coefficients <- drop(as.matrix(fit$beta)/scl$sigx)
      raw.residuals <- yy - cbind(1, xx) %*% c(a00, raw.coefficients)
      raw.rmse <- sqrt(mean(raw.residuals^2))
      raw.wt <- weight.gaussian(resi = raw.residuals, 
                                ind = indexbest, 
                                del = del)$we
      sclw <- prepara(x = xx, 
                      y = yy, 
                      family = family, 
                      index = which(raw.wt == 1), 
                      robu = 0)
      xss <- sclw$xnor
      yss <- sclw$ycen
      
      if ((missing(lambdaw))) {
        reweighted_cv <- cv.glmnet(x = xss[which(raw.wt == 1), ],
                                   y = yss[which(raw.wt == 1)], 
                                   family = family, 
                                   nfolds = 5, 
                                   alpha = alphabest, 
                                   standardize = FALSE, 
                                   intercept = FALSE, 
                                   type.measure = "mse") # NEW: REMOVED $lambda.min here because we extract it later anyways
      }
      else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw <- lambdaw
      }
      else if (!missing(lambdaw) & length(lambdaw) > 1) {
        reweighted_cv <- cv.glmnet(x = xss[which(raw.wt == 1), ],
                                   y = yss[which(raw.wt == 1)], 
                                   family = family, 
                                   lambda = lambdaw, 
                                   nfolds = 5, 
                                   alpha = alphabest, 
                                   standardize = FALSE, 
                                   intercept = FALSE, 
                                   type.measure = "mse")
      }
      
      # NEW: Choosing lambda
      if (type_lambdaw == "min") {
        lambdaw <- lambdaw_fit$lambda.min 
      } else if (type_lambdaw == "1se") {
        lambdaw <- lambdaw_fit$lambda.1se
      }
      
      # Fitting using optimal lambdaw
      fitw <- glmnet(x = xss[which(raw.wt == 1), ], 
                     y = yss[which(raw.wt == 1)], 
                     family = "gaussian", 
                     alpha = alphabest, 
                     lambda = lambdaw, 
                     standardize = FALSE, 
                     intercept = FALSE)
      a0 <- if (intercept == FALSE) 
        0
      else drop(sclw$muy + fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% (sclw$mux/sclw$sigx))
      
      coefficients <- drop(as.matrix(fitw$beta)/sclw$sigx)
      reweighted.residuals <- yy - cbind(1, xx) %*% c(a0, coefficients)
      reweighted.rmse <- sqrt(mean(reweighted.residuals^2))
      wgt <- weight.gaussian(resi = reweighted.residuals, 
                             ind = raw.wt == 1, 
                             del = del)$we
    } # End of Case: Gaussian
    # End of isTRUE(scal)
  
    
      
  #----------------------------------SCALING FALSE --------------------------------#
  # SCALING FALSE
  } else { # Think (!) this amounts to having isFALSE(scal):
    # Case: No scaling (unlikely)
    fit <- glmnet(x = x[indexbest, ], 
                  y = y[indexbest, ], 
                  family = family, 
                  alpha = alphabest, 
                  lambda = lambdabest, 
                  standardize = FALSE, 
                  intercept = FALSE)
    if (family == "binomial") {
      a00 <- if (intercept == FALSE) 
        0
      else drop(fit$a0 - as.vector(as.matrix(fit$beta)) %*% (sc$mux/sc$sigx))
      raw.coefficients <- drop(as.matrix(fit$beta)/sc$sigx)
      raw.residuals <- -(y * x %*% as.matrix(fit$beta)) + log(1 + exp(x %*% as.matrix(fit$beta)))
      raw.wt <- weight.binomial(x = xx, 
                                y = yy, 
                                beta = c(a00, raw.coefficients), 
                                intercept = intercept, 
                                del = del)
      # This is probably the default:
      if (missing(lambdaw)) {
        #reweighted_cv <- cv.glmnet(x = x[which(raw.wt == 1), ], 
        #                          y = y[which(raw.wt == 1)], 
        #                          family = family, 
        #                          nfolds = 10, 
        #                          alpha = alphabest, 
        #                          standardize = FALSE, 
        #                          intercept = FALSE, 
        #                          type.measure = "deviance") #$lambda.min # removed check check NEW
        reweighted_cv <- cva.glmnet(x = x[which(raw.wt == 1), ], 
                                    y = y[which(raw.wt == 1)], 
                                    family = family, 
                                    nfolds = 10, 
                                    #alpha = alphabest, # TODO (Check)
                                    standardize = FALSE, 
                                    intercept = FALSE)
        alphaw <- reweighted_cv$alpha[which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm)))]
        print(paste0("The reweighted alpha is: ", alphaw))
        
        
        
        
      }
      else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw <- lambdaw
      }
      # In case of multiple user supplied lambdaws
      else if (!missing(lambdaw) & length(lambdaw) > 1) {
        # NEW: saved model first, then extract the lambda we want
        reweighted_cv <- cv.glmnet(x = x[which(raw.wt == 1), ], 
                                   y = y[which(raw.wt == 1)], 
                                   family = family,
                                   lambda = lambdaw, 
                                   nfolds = 10, 
                                   alpha = alphabest, 
                                   standardize = FALSE, 
                                   intercept = FALSE, 
                                   type.measure = "deviance") # NEW: changed type.measure = "mse" to type.meaure = "deviance" 
        # NEW: Plotting the cross-validation curves
        
        
        # NEW: removed cv.glmnet(...)$lambda.min and given choice!
        if(type_lambdaw == "min") {
          lambdaw <- reweighted_fit$lambda.min
        } else if(type_lambdaw == "1se") {
          lambdaw <- reweighted_fit$lambda.1se
        }
      }
      
      # Fitting using optimal lambdaw
      fitw <- glmnet(x = x[which(raw.wt == 1), ], 
                     y = y[which(raw.wt == 1)], 
                     family = family, 
                     alpha = alphabest, 
                     lambda = lambdaw, 
                     standardize = FALSE, 
                     intercept = FALSE)
      a0 <- if(intercept == FALSE) 
        0
      else drop(fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% (sc$mux/sc$sigx))
      coefficients <- drop(as.matrix(fitw$beta)/sc$sigx)
      wgt <- weight.binomial(x = xx, 
                             y = yy, 
                             beta = c(a0, coefficients), 
                             intercept = intercept, 
                             del = del)
      reweighted.residuals <- -(yy * cbind(1, xx) %*% c(a0, coefficients)) + log(1 + exp(cbind(1, xx) %*% c(a0, coefficients)))
    } else if (family == "gaussian") { # Case: Gaussian (and no scaling)
      a00 <- if (intercept == FALSE) 
        0
      else drop(sc$muy + fit$a0 - as.vector(as.matrix(fit$beta)) %*% (sc$mux/sc$sigx))
      raw.coefficients <- drop(as.matrix(fit$beta)/sc$sigx)
      raw.residuals <- yy - cbind(1, xx) %*% c(a00, raw.coefficients)
      raw.rmse <- sqrt(mean(raw.residuals^2))
      raw.wt <- weight.gaussian(resi = raw.residuals, 
                                ind = indexbest, 
                                del = del)$we
      
      if (missing(lambdaw)) {
        reweighted_cv <- cv.glmnet(x = x[which(raw.wt == 1), ], 
                                   y = y[which(raw.wt == 1)], 
                                   family = family, 
                                   nfolds = 5, 
                                   alpha = alphabest, 
                                   standardize = FALSE, 
                                   intercept = FALSE, 
                                   type.measure = "mse") # $lambda.min # NEW REMOVED
        
      }
      else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw <- lambdaw
      }
      else if (!missing(lambdaw) & length(lambdaw) > 1) {
        reweighted_cv <- cv.glmnet(x = x[which(raw.wt == 1), ], 
                                   y = y[which(raw.wt == 1)], 
                                   family = family, 
                                   lambda = lambdaw, 
                                   nfolds = 5, 
                                   alpha = alphabest, 
                                   standardize = FALSE, 
                                   intercept = FALSE, 
                                   type.measure = "mse")
      }
      
      # NEW: lambda choosing part
      if (type_lambdaw == "min") {
        #lambdaw <- reweighted_cv$lambda.min
        lambdaw <- reweighted_cv$lambda[which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm)))]  # TODO (TEMP) Correct, this is to see if it works
        
        
      } else if (type_lambdaw == "1se") {
        #lambdaw <- reweighted_cv$lambda.1se
        lambdaw <- reweighted_cv$lambda[which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm)))]  # TODO (TEMP) correct this, this is to see if it works
        
      }
      
      # Fitting using optimal lambda
      fitw <- glmnet(x = x[which(raw.wt == 1), ], 
                     y = y[which(raw.wt == 1)], 
                     family = family, 
                     alpha = alphabest, 
                     lambda = lambdaw, 
                     standardize = FALSE, 
                     intercept = FALSE)
      a0 <- if(intercept == FALSE) 
        0
      else drop(sc$muy + fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% (sc$mux/sc$sigx))
      coefficients <- drop(as.matrix(fitw$beta)/sc$sigx)
      reweighted.residuals <- yy - cbind(1, xx) %*% c(a0, coefficients)
      reweighted.rmse <- sqrt(mean(reweighted.residuals^2))
      wgt <- weight.gaussian(resi = reweighted.residuals, 
                             ind = raw.wt == 1, 
                             del = del)$we
    }
  } # End scal == FALSE
  
  # Plotting CV Plot 
  plot(reweighted_cv)
  print(paste0("The optimal reweighted lambda is: ", lambdaw)) # Printing chosen lambdaw
  print(paste0("The optimal alpha used in the reweighting step is: ", alphabest)) # Printing chosen alpha
  # Note: I put this here at the end because otherwise I would need to repeat the statement multiple times
  
  ## PREPARING OUTPUT
  # Counting number of nonzero coefficients
  num.nonzerocoef <- sum(coefficients != 0)
  
  # Adding intercept
  intercept <- isTRUE(intercept)
  if (intercept) { 
    xx <- addIntercept(x = xx)
    coefficients <- c(a0, coefficients)
    raw.coefficients <- c(a00, raw.coefficients)
  } else {
    coefficients <- coefficients
    raw.coefficients <- raw.coefficients
  }
  if (family == "binomial") {
    u <- xx %*% raw.coefficients # XBeta?
    raw.fitted.values <- if (type == "class") {
      ifelse(test = u <= 0.5, yes = 0, no = 1)
    } else if (type == "response"){
      1/(1 + exp(-u))
    }
    uu <- xx %*% coefficients
    fitted.values <- if (type == "class") {
      ifelse(test = uu <= 0.5, yes = 0, no = 1)
    } else if (type == "response") {
      1/(1 + exp(-uu))
    }
  } else if (family == "gaussian") {
    raw.fitted.values <- xx %*% raw.coefficients
    fitted.values <- xx %*% coefficients
  }
  
  # TO DO: CHECK THIS
  if (family == "binomial") {
    objective <- h * (mean((-yy[indexbest] * (xx[indexbest, ] %*% coefficients)) + log(1 + exp(xx[indexbest, ] %*% coefficients))) + lambdabest * sum(1/2 * (1 - alphabest) * coefficients^2 + alphabest * abs(coefficients)))
  } else if (family == "gaussian") {
    objective <- h * ((1/2) * mean((yy[indexbest] - xx[indexbest, ] %*% coefficients)^2) + lambdabest * sum(1/2 * (1 - alphabest) * coefficients^2 + alphabest * abs(coefficients)))
  }
  if (intercept) {
    coefficients <- coefficients[-1] # Removing first element
    raw.coefficients <- raw.coefficients[-1]  # Removing first element
  } else {
    coefficients <- coefficients
    raw.coefficients <- raw.coefficients
  }
  inputs <- list(xx = xx, 
                 yy = yy, 
                 family = family, 
                 alphas = alphas, 
                 lambdas = lambdas, 
                 lambdaw = lambdaw, 
                 hsize = hsize, 
                 h = h, 
                 nsamp = nsamp, 
                 s1 = s1, 
                 nCsteps = nCsteps, 
                 nfold = nfold, 
                 intercept = intercept, 
                 repl = repl, 
                 para = para, 
                 ncores = ncores, 
                 del = del, 
                 scal = scal)
  
  ## OUTPUTS
  # Case: binomial
  if(family == "binomial"){
    output <- list(objective = objective, 
                   best = sort(indexbest), 
                   raw.wt = raw.wt, 
                   wt = wgt, 
                   a00 = a00, 
                   raw.coefficients = raw.coefficients, 
                   a0 = a0, 
                   coefficients = coefficients, 
                   alpha = alphabest, 
                   lambda = lambdabest, 
                   lambdaw = lambdaw, 
                   num.nonzerocoef = num.nonzerocoef, 
                   h = h, 
                   raw.residuals = drop(raw.residuals), 
                   residuals = drop(reweighted.residuals), 
                   fitted.values = drop(fitted.values), 
                   raw.fitted.values = drop(raw.fitted.values), 
                   classnames = classnames, 
                   classsize = ntab, 
                   inputs = inputs, 
                   indexall = indexall, 
                   call = sys.call(),
                   alphas = alphas,  # NEW: Added the alphas that were used
                   lambdas = lambdas)  # NEW: Added the lambdas that were used
    
  } else if (family == "gaussian"){
    output <- list(objective = objective, 
                   best = sort(indexbest), 
                   raw.wt = raw.wt, 
                   wt = wgt, 
                   a00 = a00, 
                   raw.coefficients = raw.coefficients, 
                   a0 = a0, 
                   coefficients = coefficients, 
                   alpha = alphabest, 
                   lambda = lambdabest, 
                   lambdaw = lambdaw, 
                   num.nonzerocoef = num.nonzerocoef, 
                   h = h, 
                   raw.residuals = drop(raw.residuals), 
                   residuals = drop(reweighted.residuals), 
                   fitted.values = drop(fitted.values), 
                   raw.fitted.values = drop(raw.fitted.values), 
                   raw.rmse = raw.rmse, 
                   rmse = reweighted.rmse, 
                   inputs = inputs, 
                   indexall = indexall, 
                   call = sys.call(),
                   alphas = alphas,  # NEW: Added the alphas that were used
                   lambdas = lambdas)  # NEW: Added the lambdas that were used
  }
  class(output) <- "enetLTS"
  output$call <- matchedCall
  output
}