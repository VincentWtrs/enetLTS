ic_penalty <- function(type, model, X, alpha, EBIC_sigma = 0.25, s = NULL, intercept = NULL){
  # This function calculates the PENALTY factor associated with information criteria (which then needs to be ADDED to loglik(= minus loss))
  
  ## INPUTS:
  # glmnet_in: a model fitted through glmnet()
  # type: the sort of information criterion wanted
  # X: the data
  # alpha: because for some stupid reason one cannot extract it from the glmnet function.
  # intercept: if intercept was used to fit the logit-glmnet (family = "binomial") model: TRUE.
  
  ## OUPUTS:
  # Numeric value denoting the effective degrees of freedom
  
  ## Error Catching (why commented out?)
  # Checking if model is glmnet
  #if(class(glmnet_in) != c("lognet", "glmnet")){
  #  stop("The input for ic_penalty() should be a glmnet model")
  #}
  
  if (length(model$lambda) > 1) {
    if(length(s) > 1){
      stop("The ic_penalty function cannot handle multiple lambdas without supplying a single value for lambda through the 's' argument!")
    }
  }
  
  # If the model was a full path of regularization, but a single s (lambda value) was supplied
  if(length(model$lambda) > 1) {
    if(length(s) == 1) {
      lambda <- s
    }
  }
  
  # If no s given, extract s from the model object
  if(is.null(s)) {
    lambda = model$lambda
  }
  
  # Checking if type is supported
  if(type != "AIC"  &
     type != "AIC_C" &
     type != "BIC" &
     type != "EBIC" &
     type != "GIC" & 
     type != "BIC_WLL" &
     type != "BIC_FT" &
     type != "HBIC" &
     type != "BIC_HD" &
     type != "EBIC2" &
     type != "ERIC"){
    stop("This type of information criterion not supported (or check for typos)")
  }
  
  # Extracting parameters
  nobs <- model$nobs
  
  ## Intercept can be given to function but we will perform a check if thats consistent with what the model has been fitted on
  # If an intercept is given with the function call
  if(!is.null(intercept)){
    if(!intercept) { # If intercept FALSE
      if(any(model$a0 != 0)){ # If there is one not equal to 0
        stop("There is a discrepancy between the function call to ic_penalty and the given model, since the model has a fitted intercept, while the function was called with intercept=FALSE")
      }
    }
    if(intercept) { # If intercept TRUE
      if(all(model$a0 == 0)) {
        stop("There is a discrepancy between the function call to ic_penalty and the given model, since the model appears to not have any fitted intercept but the function call states intercept=TRUE")
      }
    }
  }
  
  # Checking if intercept is in the model
  if(all(model$a0 == 0)){
    intercept <- FALSE
  } else {
    intercept <- TRUE
  }
  
  if(intercept == FALSE){
    p <- length(coefficients(model, s = s)) - 1 # Because intercept is always returned when using coefficients() (hence - 1)
  }
  # Calculating the effective degrees of freedom
  df <- logit_df(model = model,
                  X = X,
                  alpha = alpha,
                  lambda = lambda) # NEW
  
  # Amount of nonzeros
  nonzeros <- sum(coefficients(model, s = lambda) != 0)
  
  ## Calculating the information criteria penalty
  # Akaike Information Criterion (AIC)
  if(type == "AIC"){
    penalty <- (2 * df)/nobs
  }
  
  # Corrected AIC (AIC_C)
  if(type == "AIC_C"){
    penalty <- (2 * df * (nobs/(nobs - df - 1)))/nobs
  }
  
  # Bayesian Information Criterion (BIC)
  if(type == "BIC"){
    penalty <- (df * log(nobs))/nobs
  }
  
  # BIC-type: Wang, Li & Leng (2009)
  if(type == "BIC_WLL"){
    penalty <- (df * log(nobs) * log(log(p)))/nobs
  }
  
  # BIC-type: Fang & Tang (2012) -> GIC
  if(type == "BIC_FT"){ 
    penalty <- (df * log(log(n)) * log(p))/nobs
  }
  
  # HBIC: High-dimensional BIC
  if(type == "HBIC"){
    sigma <- 1.5 # Empirically well-performing ]1, 2]
    if(sigma > 1){
      stop("Sigma is required to be > 1 following Wang & Zhu (2011)")
    }
    penalty <- (2 * sigma * df * log(p))/nobs
  } # I don't see how this comes
  
  # BIC for High Dimensions (Gao, ...)
  if(type == "BIC_HD"){
    c <- 1 # Common choices: 1 or 2
    penalty <- (c * log(p) * df)/nobs
  }
  
  # Extended Bayesian Information Criterion (EBIC)
  if(type == "EBIC"){
    sigma <- EBIC_sigma # Some default value! (e.g. 0.25)
    #penalty <- (coefs_nonzero * log(nobs) + 2 * coefs_nonzero * sigma * log(p))/nobs # THEORY
    penalty <- (df * log(nobs) + 2 * df * sigma * log(p))/nobs
  }
  
  # EBIC2: Alternative EBIC - Chen & Chen (2008)
  if(type == "EBIC2"){
    sigma <- EBIC_sigma # 0.25 good default
    penalty <- (df * log(nobs) + 2 * sigma * log(choose(p, nonzeros)))/nobs
  }
  
  # Generalized Information Criterion (GIC)
  if(type == "GIC"){
    a_n <- log(log(nobs)) * log(p) # Also possible: log(p)
    penalty <- (a_n * df)/nobs
  }
  
  # Extended Regularized Information Criterion (ERIC)
  if(type == "ERIC"){
    nu <- 0.5 # Often suggested values: 0.5 or 1
    penalty <- (2 * nu * nonzeros * log(nobs/lambda))/nobs
    # Alternatively: penalty <-  (2 * nu * df * log(nobs/lambda))/nobs
  }
  
  return(penalty)
}
