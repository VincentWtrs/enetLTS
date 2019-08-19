cv.enetLTS <- function(index = NULL, xx, yy, family, h, alphas, lambdas, nfold, 
                       repl, ncores, plot = TRUE, ic_type = NULL){
  
  MNLL <- TMNLL <- RTMSPE <- RMSPE <- NULL
  n <- nrow(xx)
  p <- ncol(xx)
  
  ## Input parameter checks
  # alphas sequence
  wh <- (alphas < 0 | alphas > 1) # This is already done on the higher level function...
  if (sum(wh) > 0) 
    stop("alphas can take the values only between 0 and 1")
  if (missing(alphas)) 
    stop("provide an alphas sequence")
  
  # lambdas sequence
  if (missing(lambdas)) 
    stop("provide an lambdas sequence")
  
  evalCrit <- matrix(NA, nrow = length(lambdas), ncol = length(alphas))
  dimnames(evalCrit) <- list(paste("lambdas", lambdas), paste("alpha", alphas))
  combis_ind <- expand.grid(1:length(lambdas), 1:length(alphas))
  indcombi <- 1:nrow(combis_ind)
  
  ## REMOVED: DEFINITION OF calc_evalCrit(), THIS HAS BEEN SPLIT OFF TO SEPARATE FILE
  
  # Running for all combinations of alpha, lambda
  ### NEW: added argument family = family and ic_type to be passed on as well
  if(is.null(ic_type) | length(ic_type) == 1) {
    temp_result <- mclapply(1:nrow(combis_ind), FUN = calc_evalCrit, 
                            combis_ind = combis_ind, 
                            alphas = alphas, 
                            lambdas = lambdas, 
                            index = index, 
                            xx = xx, 
                            yy = yy, 
                            nfold = nfold, 
                            repl = repl, 
                            family = family, # NEW 
                            ic_type = ic_type, # NEW
                            mc.cores = ncores, 
                            mc.allow.recursive = FALSE)
    print("I am printing temp_result in cv.enetLTS")
    print(temp_result)
    
    # Restructuring output
    temp_result2 <- matrix(unlist(temp_result), 
                           ncol = repl + 2, # + 2
                           byrow = TRUE)
    print("printing temp_result2 in cv.enetLTS")
    print(temp_result2)
    
    
    for (k in 1:nrow(temp_result2)) {
      i <- temp_result2[k, 1]
      j <- temp_result2[k, 2]
      evalCrit[i, j] <- mean(temp_result2[k, 3:(repl + 2)])
    }
    
    print("I am printing evalCrit in cv.enetLTS:")
    print(evalCrit)
    
    optind <- which(evalCrit == min(evalCrit, na.rm = TRUE), arr.ind = TRUE)[1, ]
    minevalCrit <- evalCrit[optind[1], optind[2]]
    indexbest <- index[, optind[1], optind[2]]
    #alphas <- round(alphas, 10) # NEW: REMOVED ROUNDING
    alpha <- alphas[optind[2]] 
    #lambdas <- round(lambdas, 10) # NEW: less rounding # NEW: REMOVED ROUNDING
    lambda <- lambdas[optind[1]] # NEW: less rounding
    
    ## PLOTTING
    # NEW: was fully defined inside this function but split off to new function
    if(isTRUE(plot)){
      tune_plot(alphas = alphas,
                alpha = alpha,
                lambdas = lambdas,
                lambda = lambda,
                evalCrit = evalCrit,
                index = index,
                family = family)
    }
  } else if (length(ic_type) > 1) {
    temp_result <- mclapply(1:nrow(combis_ind), FUN = calc_evalCrit2, 
                            combis_ind = combis_ind, 
                            alphas = alphas, 
                            lambdas = lambdas, 
                            index = index, 
                            xx = xx, 
                            yy = yy, 
                            nfold = nfold, 
                            repl = repl, 
                            family = family, # NEW 
                            ic_type = ic_type, # NEW
                            mc.cores = ncores, 
                            mc.allow.recursive = FALSE)
    print("I am printing temp_result in cv.enetLTS")
    print(temp_result)
    
    # NEW UNPACKING
    for(i in 1:length(ic_type)) {
      temp_result2[[i]] <- matrix(unlist(temp_result[[i]]),
                                  ncol = repl + 2,
                                  byrow = TRUE)
    }
    
    # Restructuring output
    #temp_result2 <- matrix(unlist(temp_result), 
    #                       ncol = repl + 2, # + 2
    #                       byrow = TRUE)
    print("printing temp_result2 in cv.enetLTS")
    print(temp_result2)
    
    for(m in 1:length(ic_type)){
      for (k in 1:nrow(temp_result2)) {
        i <- temp_result2[k, 1]
        j <- temp_result2[k, 2]
        evalCrit[m, i, j] <- mean(temp_result2[k, 3:(repl + 2)]) # Array
      }
    }

    
    print("I am printing evalCrit in cv.enetLTS:")
    print(evalCrit)
    
    optind <- which(evalCrit == min(evalCrit, na.rm = TRUE), arr.ind = TRUE)[1, ]
    minevalCrit <- evalCrit[optind[1], optind[2]]
    indexbest <- index[, optind[1], optind[2]]
    #alphas <- round(alphas, 10) # NEW: REMOVED ROUNDING
    alpha <- alphas[optind[2]] 
    #lambdas <- round(lambdas, 10) # NEW: less rounding # NEW: REMOVED ROUNDING
    lambda <- lambdas[optind[1]] # NEW: less rounding
    
    ## PLOTTING
    # NEW: was fully defined inside this function but split off to new function
    if(isTRUE(plot)){
      tune_plot(alphas = alphas,
                alpha = alpha,
                lambdas = lambdas,
                lambda = lambda,
                evalCrit = evalCrit,
                index = index,
                family = family)
    }
  }
  
  # OUTPUT
  return(list(evalCrit = evalCrit, 
              minevalCrit = minevalCrit, 
              indexbest = indexbest, 
              lambdaopt = lambda, 
              alphaopt = alpha))
}