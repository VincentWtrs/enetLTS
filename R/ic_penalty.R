ic_penalty <- function(type, model, X, alpha, intercept, EBIC_sigma = 0.25){
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
  
  # Amount of parameters (p of the design matrix basically)
  if(intercept == TRUE){
    stop("IC with intercept currently not implemented yet")
  }
  if(intercept == FALSE){
    p <- length(coefficients(model)) - 1 # Because intercept is always returned (hence - 1)
  }
  
  # Calculating Degrees of Freedom (Effective DF!)
  df <- logit_df(model = model,
                 X = X,
                 alpha = alpha,
                 intercept = intercept)
  
  # Amount of nonzeros
  nonzeros <- model$df
  
  ## Calculating the information criteria penalty
  # Akaike Information Criterion (AIC)
  if(type == "AIC"){
    penalty <- (2 * df)/nobs
  }
  
  # Corrected AIC (AIC_C)
  if(type == "AIC_C"){
    penalty <- (2 * df * (nobs/(n - df - 1)))/nobs
  }
  
  # Bayesian Information Criterion (BIC)
  if(type == "BIC"){
    penalty <- (df * log(nobs))/nobs
  }
  
  # BIC-type: Wang, Li & Leng (2009)
  if(type == "BIC_WLL"){
    penalty <- (df * log(n) * log(log(p)))/nobs
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
