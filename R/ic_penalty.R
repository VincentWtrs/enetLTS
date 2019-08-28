ic_penalty <- function(type, model, X, alpha, naive = FALSE, HBIC_gamma = NULL, BIC_HD_c = NULL, EBIC_gamma = NULL, 
                       ERIC_nu = NULL, s = NULL, intercept = NULL, brute_force = FALSE) {
  # This function calculates the PENALTY factor associated with information criteria (which then needs to be ADDED to loglik(= minus loss))
  
  ## INPUTS:
  # type: The specific information criterion requested
  # model: A glmnet model, either a model with a single lambda or a regularization path (multiple lambdas)
  # X: The design matrix without the intercept!
  # alpha: The Elastic Net mixing hyperparameter
  # naive: If naive calculations using the amount of estimated parameters is to be used for effective df (Default: FALSE)
  # HBIC_gamma: Tuning parameter of the HBIC IC
  # BIC_HD_c: Tuning parameter of the BIC_HD Information Criterion
  # EBIC_gamma: Tuning parameter of the Extended BIC
  # ERIC_nu: Tuning parameter of the Extended Regularized IC
  # s: Single value for regularization parameter lambda in case where argument 'model' contained a regularization path
  # intercept: Optional parameter for clarity: shows if intercept has been fitted or not and will adjust edfs appropriately. Will check against provided model
  
  ## OUPUTS:
  # Numeric value denoting the effective degrees of freedom
  
  ## TODO 
  # 1. Make the naming of lambda / s clearer in both logit_df as well as ic_penalty
  # 2. Include possibility to estimate nuisance parameters (e.g. add amount of nuisance parameters as integer argument to function, ...)
  # 3. Introduce the detect_path() and detect_intercept() functions that are in the logit_df function.
  
  ## Error Catching
  # Checking if model is glmnet
  if (class(model)[1] != "lognet" & class(model)[2] != "glmnet") {
    stop("The input for ic_penalty() should be a glmnet model")
  }
  
  
  # Checking for s is regularization path model was supplied
  if (length(model$lambda) > 1) {
    if (length(s) > 1) {
      stop("The ic_penalty function cannot handle multiple lambdas without supplying a single value for lambda through the 's' argument!")
    }
  }
  
  # Intercept can be given to function but we will perform a check if thats consistent with what the model has been fitted on
  if (!is.null(intercept)) {  # If an intercept is given with the function call
    if (isFALSE(intercept)) { # If intercept FALSE
      if (any(model$a0 != 0)) {
        # In this case the user requests no intercept, but intercepts have been fitted! (Using any() to check in case of regularization path)
        stop("There is a discrepancy between the function call to ic_penalty and the given model, since the model has a fitted intercept, while the function was called with intercept=FALSE")
      }
    } else if (isTRUE(intercept)) {
      if (all(model$a0 == 0)) {
        # In this case the uer says there is an intercept but none have been fitted (Using any() to check in case of regularization path)
        stop("There is a discrepancy between the function call to ic_penalty and the given model, since the model appears to not have any fitted intercept but the function call states intercept=TRUE")
      }
    }
  }
  
  # Checking if IC type is supported
  if(type != "AIC"  &
     type != "AIC_C" &
     type != "BIC" &
     type != "EBIC" &
     type != "GIC" & 
     type != "BIC_WLL" &
     type != "GIC_FT" &
     type != "HBIC" &
     type != "BIC_HD" &
     type != "EBIC2" &
     type != "ERIC" &
     type != "RIC") {
    stop("This type of information criterion not supported (or check for typos)")
  }
  
  ## Extracting parameters from model
  # Extracting lambda
  if(is.null(s)) { # If no s given: get it from model itself
    lambda = model$lambda
  } else if (length(model$lambda) > 1) { # If s is given, pass it on
    lambda <- s 
    # Note: no additional check on length(s) needs to be done because would have thrown stop() above
  }
  
  # Extracting amoutn of observations
  nobs <- model$nobs # Using nobs instead of n because here the amount of observations will often be |H| = h
  
  # Extracting intercept parameter
  if (all(model$a0 == 0)) {
    intercept <- FALSE
  } else if (!all(model$a0 == 0)) {
    intercept <- TRUE
  }
  
  # Extracting total available predictor dimensionality (p)
  if (intercept == FALSE) {
    p <- length(coefficients(model, s = lambda)) - 1 # Because intercept is always returned when using coefficients() (hence - 1)
  } else if (isTRUE(intercept)) {
    p <- length(coefficients(model, s = lambda)) # In this case coefficients() returns the correct amount
  }
  # TODO: Maybe use p always and then make later different flows to use p+1 instead of p, now the definition of p changes depending on intercept
  
  # Calcualting amount of nonzero predictors
  nonzeros <- sum(coefficients(model, s = lambda) != 0)
  # NOTE: No different path for intercept TRUE/FALSE needed: if no intercept fitted: it will be 0 anyway, if there is one fitted it will be in the p
  
  # Calculating the effective degrees of freedom (df)
  if(isFALSE(naive)) {
    df <- logit_df(model = model,
                   X = X,
                   alpha = alpha,
                   lambda = lambda, # NEW
                   brute_force = brute_force) # NEW
    # TODO: CHECK DOES IT HANDLE INTERCEPT WELL?
  } else if (isTRUE(naive)) {
    df <- nonzeros # If naive calculation is asked just give the amount of nonzeros
  }
  
  ## Calculating the information criteria penalty
  # Akaike Information Criterion (AIC)
  if (type == "AIC") {
    penalty <- (2 * df)/nobs
  }
  
  # Corrected AIC (AIC_C)
  if (type == "AIC_C") {
    penalty <- (2 * df * (nobs/(nobs - df - 1)))/nobs
  }
  
  # Bayesian Information Criterion (BIC)
  if (type == "BIC") {
    penalty <- (df * log(nobs))/nobs
  }
  
  # BIC-type: Wang, Li & Leng (2009)
  if (type == "BIC_WLL") {
    penalty <- (df * log(nobs) * log(log(p)))/nobs
  }
  
  # GIC-type: Fang & Tang (2012) -> GIC
  if (type == "GIC_FT") { 
    penalty <- (df * (log(log(nobs)) * log(p)))/nobs
  }
  
  # HBIC: High-dimensional BIC
  if (type == "HBIC") {
    if (is.null(HBIC_gamma)) {
      HBIC_gamma <- 1.5  # Empirically well-performing ]1, 2]
    }
    if (HBIC_gamma < 1) {
      warning("HBIC_gamma ideally set above 1 to become selection consistent")
    }
    penalty <- (2 * HBIC_gamma * df * log(p))/nobs
  } 
  
  # BIC for High Dimensions (Gao & Carroll (2017))
  if(type == "BIC_HD") {
    if (is.null(BIC_HD_c)) {
      BIC_HD_c <- 1 # Common choices: 1 or 2
    }
    penalty <- (BIC_HD_c * log(p) * df)/nobs
  }
  
  # Extended Bayesian Information Criterion (EBIC)
  if (type == "EBIC") {
    if (is.null(EBIC_gamma)) {
      EBIC_gamma <- 0.25
    }
    if(EBIC_gamma > 1 | EBIC_gamma < 0) {
      stop("This EBIC_gamma parameter is not allowed, choose a value in [0, 1]")
    }
    #penalty <- (coefs_nonzero * log(nobs) + 2 * coefs_nonzero * sigma * log(p))/nobs # THEORY
    penalty <- (df * log(nobs) + 2 * df * EBIC_gamma * log(p))/nobs
  }
  
  # EBIC2: Alternative EBIC - Chen & Chen (2008)
  if (type == "EBIC2") {
    if (is.null(EBIC_gamma)) {
      EBIC_gamma <- 0.25 # 0.25 good default
    }
    penalty <- (df * log(nobs) + 2 * EBIC_gamma * log(choose(p, nonzeros)))/nobs
  }
  
  # Generalized Information Criterion (GIC)
  if(type == "GIC"){ 
    a_n <- log(log(nobs)) * log(p) # Also possible: log(p)
    penalty <- (a_n * df)/nobs
  }
  
  # Extended Regularized Information Criterion (ERIC)
  if (type == "ERIC") {
    if (is.null(ERIC_nu)) {
      ERIC_nu <- 0.5 # Often suggested values: 0.5 or 1
    }
    penalty <- (2 * ERIC_nu * nonzeros * log(nobs/lambda))/nobs
    # Alternatively: penalty <-  (2 * nu * df * log(nobs/lambda))/nobs
  }
  
  # RIC
  if (type == "RIC") {
    penalty <- (2 * log(p) * df)/nobs # Simplified version of HBIC (HBIC_gamma == 1)
  }
  
  # OUTPUT
  return(penalty)
}
