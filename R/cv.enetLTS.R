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
  if (length(ic_type) > 1) {
    temp_result <- mclapply(1:nrow(combis_ind), FUN = calc_evalCrit2, # CALLING calc_evalCrit2 (adjusted function) 
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
    ## Handling results
    # Make an expanded grid
    temp_grid <- expand.grid("IC" = ic_type, "alpha" = alphas, "lambda" = lambdas, "loss" = NA)
    
    # Filling grid
    s <- 1 # Counter
    for(r in 1:length(temp_result)) {
      for(m in 1:length(ic_type)) {
        temp_grid[s, "lambda"] <- lambdas[temp_result[[r]]$lambda_ind]
        temp_grid[s, "alpha"] <- alphas[temp_result[[r]]$alpha_ind]
        temp_grid[s, "IC"] <- ic_type[m]
        temp_grid[s, "loss"] <- temp_result[[r]]$evalCritl[m]
        s <-  s + 1
      }
    }
    print("printing temp_grid")
    print(temp_grid)
    
    ## Best for each IC
    # Initializations
    best <- rep(NA, times = length(ic_type))
    alpha_opt <- rep(NA, times = length(ic_type))
    lambda_opt <- rep(NA, times = length(ic_type))
    loss_opt <- rep(NA, times = length(ic_type))
    output <- vector("list", length = length(ic_type))
    alpha_indx <- rep(NA, times = length(ic_type))
    lambda_indx <- rep(NA, times = length(ic_type))
    
    
    # Loop
    for(m in 1:length(ic_type)) {
      print("printing ic_Type[m]:")
      print(ic_type[m])
      grid_ic_now <- temp_grid[temp_grid$IC == ic_type[m], ]
      print("Printing grids_ic_now:")
      print(grid_ic_now)
      alpha_opt[m] <- grid_ic_now[which.min(grid_ic_now$loss), "alpha"]
      lambda_opt[m] <- grid_ic_now[which.min(grid_ic_now$loss), "lambda"]
      loss_opt[m] <- grid_ic_now[which.min(grid_ic_now$loss), "loss"]
      alpha_indx[m] <- which(alphas == alpha_opt[m])
      lambda_indx[m] <- which(lambdas == lambda_opt[m])
      
      output[[m]]$minevalCrit <- loss_opt[m]
      
      print("printing the alpha_indx & lambdaindx")
      print("for index m")
      print(alpha_indx[m])
      print(lambda_indx[m])
      
      
      output[[m]]$indexbest <- index[, lambda_indx[m], alpha_indx[m]] # NOTE THE OPPOSITE INDEXATION (WEIRD BUT CORRECT) TODO DOUBLE CHECK
      output[[m]]$alphas <- alphas
      output[[m]]$lambdas <- lambdas
      output[[m]]$alpha <- alpha_opt[m]
      output[[m]]$lambda <- lambda_opt[m]
    }
    # OUTPUT IS A LIST NOW
    print("best alpha/lambda indices")
    print(alpha_indx)
    print(lambda_indx)
    
    print("Best values of alpha/lambda")
    print(alphas[alpha_indx])
    print(lambdas[lambdas_indx])
    
    return(output)
    
    
    
  } else { # IN CASE NO MULTIPLE ICS! (CV or single IC)
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
    
    # Restructuring output
    temp_result2 <- matrix(unlist(temp_result), 
                           ncol = repl + 2, # + 2
                           byrow = TRUE)
    
    for (k in 1:nrow(temp_result2)) {
      i <- temp_result2[k, 1]
      j <- temp_result2[k, 2]
      evalCrit[i, j] <- mean(temp_result2[k, 3:(repl + 2)])
    }
    
    optind <- which(evalCrit == min(evalCrit, na.rm = TRUE), arr.ind = TRUE)[1, ]
    minevalCrit <- evalCrit[optind[1], optind[2]]
    indexbest <- index[, optind[1], optind[2]]
    alphas <- round(alphas, 10)
    alpha <- alphas[optind[2]] 
    lambdas <- round(lambdas, 10) # NEW: less rounding
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
    
    # OUTPUT
    return(list(evalCrit = evalCrit, 
                minevalCrit = minevalCrit, 
                indexbest = indexbest, 
                lambdaopt = lambda, 
                alphaopt = alpha))
  }
}