CStep <- function(x,y,family,indx,h,hsize,alpha,lambda,scal){
  ## internal function
  
  # require(glmnet)
  # source("utilities.R")
  # source("objectiveFunc.R")
  
  n <- nrow(x)
  if (scal){
    scl <- prepara(x = x,
                   y = y,
                   family = family,
                   index = indx,
                   robu = 0)
    xs <- scl$xnor
    ys <- scl$ycen
    if (family == "binomial") {
      fit <- glmnet(xs[indx,],ys[indx],family,alpha=alpha,lambda=lambda,standardize=FALSE,intercept=FALSE)
      beta <- matrix(fit$beta)
      resid <- -(ys * xs %*% beta) + log(1 + exp(xs %*% beta))
      if (all(beta == 0)){
        return(list(object = -Inf,index = indx,residu = resid, beta = beta))
        } 
      resid.sort <- sort(resid, decreasing = FALSE, index.return = TRUE) 
      h0 <- floor((length(y[y == 0]) + 1) * hsize)
      h1 <- h - h0
      index0 <- resid.sort$ix[y[resid.sort$ix] == 0][1:h0]
      index1 <- resid.sort$ix[y[resid.sort$ix] == 1][1:h1]
      indxnew <- c(index0, index1)
    } else if(family == "gaussian"){
      fit <- glmnet(xs[indx,],ys[indx],family,alpha=alpha,lambda=lambda,standardize=FALSE,intercept=FALSE)
      beta <- matrix(fit$beta)
      resid <- ys - predict(fit,xs,exact=TRUE)
      resid.sort <- sort(abs(resid),index.return=TRUE)
      indxnew <- resid.sort$ix[1:h]
    }
    obj <- Objval(x = xs,
                  y = ys,
                  family = family,
                  coef = beta,
                  ind = indxnew,
                  alpha = alpha,
                  lambda = lambda)
  
  } else {
    if (family == "binomial") {
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
        return(list(object = -Inf,index = indx,residu=resid,beta=beta))}
      resid.sort <- sort(resid, decreasing = FALSE, index.return = TRUE) 
      h0 <- floor((length(y[y == 0]) + 1) * hsize)
      h1 <- h-h0
      index0 <- resid.sort$ix[y[resid.sort$ix] == 0][1:h0]
      index1 <- resid.sort$ix[y[resid.sort$ix]==1][1:h1]
      indxnew <- c(index0,index1)
      
    } else if (family == "gaussian") {
      fit <- glmnet(x = x[indx,],
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
  return(list(object=obj,index=indxnew,residu=resid,beta=beta))
}