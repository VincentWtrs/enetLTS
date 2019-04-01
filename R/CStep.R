CStep <- function(x, y, family, indx, h, hsize, alpha, lambda, scal){
  ## internal function
  
  # require(glmnet)
  # source("utilities.R")
  # source("objectiveFunc.R")
  
  n <- nrow(x)
  
  # Case: Scaling data (nonrobustly) (DEFAULT)
  if (scal){
    scl <- prepara(x = x,
                   y = y,
                   family = family,
                   index = indx, # TO DO: Can the data be scaled nonrobustly because it is already outlier-free?
                   robu = 0)
    xs <- scl$xnor
    ys <- scl$ycen
    
    # Case: Scaling & Binomial
    if (family == "binomial") {
      # Fitting Elastic Net with given settings
      #fit <- glmnet(x = xs[indx,],
      #              y = ys[indx],
      #              family = family,
      #              alpha = alpha, # Single alpha
      #              lambda = lambda, # Single lambda
      #              standardize = FALSE,
      #              intercept = FALSE) # DUBIOUS! If changing this also change beta and resid! BUT KEEP TRACK WITH THE if(all(beta==0))
      
      # NEW: with intercept = TRUE
      fit <- glmnet(x = xs[indx, ],
                    y = ys[indx],
                    family = family,
                    alpha = alpha,
                    lambda = lambda,
                    standardize = FALSE,
                    intercept = TRUE) # NEW set to TRUE
      
      beta_with_int <- matrix(coef(fit)) # NEW: to include b0
      beta <- matrix(fit$beta) # Getting beta ($beta gets coefs WITHOUT INTERCEPT!) # We can still use it for the test later
      #resid <- -(ys * xs %*% beta_with_int) + log(1 + exp(xs %*% beta_with_int) # OLD
      resid <- -(ys * cbind(1, xs) %*% beta_with_int) + log(1 + exp(cbind(1, xs) %*% beta_with_int)) # NEW: used beta_with_int and cbind(1, Xs) to accomodate
      
      
      
      # Fallback if all beta == 0 # Stop early (?)
      if (all(beta == 0)){
        return(list(object = -Inf,index = indx, residu = resid, beta = beta))
      }
      
      # If not all equal to 0: 
      resid.sort <- sort(resid, decreasing = FALSE, index.return = TRUE) 
      h0 <- floor((length(y[y == 0]) + 1) * hsize)
      h1 <- h - h0
      index0 <- resid.sort$ix[y[resid.sort$ix] == 0][1:h0]
      index1 <- resid.sort$ix[y[resid.sort$ix] == 1][1:h1]
      indxnew <- c(index0, index1)
    } else if (family == "gaussian") {
      fit <- glmnet(x = xs[indx,],
                    y = ys[indx],
                    family = family,
                    alpha = alpha,
                    lambda = lambda,
                    standardize = FALSE,
                    intercept = FALSE)
      beta <- matrix(fit$beta)
      resid <- ys - predict(fit,xs,exact=TRUE)
      resid.sort <- sort(abs(resid),index.return=TRUE)
      indxnew <- resid.sort$ix[1:h]
    }
    
    # Calculating objective function
    obj <- Objval(x = xs,
                  y = ys,
                  family = family,
                  coef = beta,
                  ind = indxnew, # Note: it gets the index supplied!
                  alpha = alpha,
                  lambda = lambda)
    
  # Case: No scaling (Unlikely)
  } else if (isFALSE(scal)) {
    if (family == "binomial") { # TO DO: Correct this later
      fit <- glmnet(x = x[indx, ],
                    y = y[indx],
                    family = family,
                    alpha=  alpha,
                    lambda = lambda,
                    standardize = FALSE,
                    intercept = FALSE)
      beta <- matrix(fit$beta)
      resid <- -(y * x %*% beta) + log(1+exp(x %*% beta))
      if (all(beta == 0)) {
        return(list(object = -Inf,index = indx,residu =resid,beta=beta))}
      resid.sort <- sort(resid, decreasing = FALSE, index.return = TRUE) 
      h0 <- floor((length(y[y == 0]) + 1) * hsize)
      h1 <- h - h0
      index0 <- resid.sort$ix[y[resid.sort$ix] == 0][1:h0]
      index1 <- resid.sort$ix[y[resid.sort$ix] == 1][1:h1]
      indxnew <- c(index0,index1)
      
    } else if (family == "gaussian") {
      fit <- glmnet(x = x[indx, ],
                    y = y[indx],
                    family = family,
                    alpha = alpha,
                    lambda = lambda,
                    standardize = FALSE,
                    intercept = FALSE)
      beta <- matrix(fit$beta)
      resid <- y - predict(fit,x, exact = TRUE)
      resid.sort <- sort(abs(resid), index.return = TRUE)
      indxnew <- resid.sort$ix[1:h]
    }
    
    obj <- Objval(x = x,
                  y = y,
                  family = family,
                  coef = beta,
                  ind = indxnew,
                  alpha = alpha,
                  lambda = lambda)
  }
  
  # OUTPUT
  return(list(object=obj,index=indxnew,residu=resid,beta=beta))
}