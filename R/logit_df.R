logit_df <- function(model, X, alpha = NULL, s = NULL){
  ### logit_df() FUNCTION: calculating the effective degrees of freedom in a potentially penalized logistic regression using hat matrix trace
  
  ## INPUTS
  # model: A glmnet() model
  # X: the data matrix WITHOUT INTERCEPT
  # alpha: Alpha value used to fit the model (single value)
  # intercept: If intercept was used to fit the logit-glmnet model: set to TRUE
  
  ## OUTPUT
  # df: the effective degrees of freedom for the logit model
  
  ## NOTES
  # 1. lamba value can be extracted from glmnet model, alpha cannot in most circumstances
  # 2. The case for including an intercept is not so clean as it seems as it would require different
  #    inner workings for the predict() function as it does not work with matrices with intercepts... might do this later
  # 3. IDEA to work with the intercept is to see in coefficients(model) is it's 0, if yes: no intercept fitted, becuase it's always unpenalized!
  # 4. I was afraid this might break down at the matrix inversion when |Active set| > n, but seems to be okay!
  
  detect_intercept <- function(model) {
    if(any(model$a0 !=  0)){
      intercept <- TRUE
    } else {
      intercept <- FALSE
    }
    # OUTPUT
    return(intercept)
  }
  
  detect_regularization_path <- function(model) {
    if(length(model$lambda) > 1){
      path <- TRUE
    } else {
      path <- FALSE
    }
  }
  
  # Checking if the model was fitted with an intercept
  intercept <- detect_intercept(model)
  
  # Checking if the model was fitted over a regularization path
  path <- detect_regularization_path(model)
  
  # Checking, if there was a regularization path fitted that an additional single labmda is supplied and then passing on that lambda if so
  if(path) {
    if(is.null(s)){
      stop("logit_df was supplied a glmnet model fitted over a regularization path but no specific lambda was given to calculate the degrees of freedom")
    } else {
      # Setting the supplied lambda as the lambda to further use
      if(length(s) == 1){ # Additional check, we can only have a single lambda here
        lambda <- s
      }
    }
  }
  
  # Trying to retrieve alpha if missing in function call
  if(missing(alpha)){
    if(is.numeric(model$call$alpha)){
      alpha <- model$call$alpha # Trying to extract from the function call if possible
    } else { # But if it's another variable, it will return e.g. "myalpha" and will not be of any use
      stop("No alpha value given for logit_df(), also no value could be succesfully extracted from the model object, supply a valid alpha: [0, 1]")
    }
  }
  
  ## Construction of the design matrix
  # If there is an intercept, add a first column of 1
  if(intercept) {
    X <- cbind(1, X)
    colnames(X)[1] <- "X0"
  } # If there is no intercept, the X matrix is good as is
  
  ## Extracting fitted beta coefficients
  # If intercept, get the beta0 in there as well!
  if(intercept) {
    beta <- as.matrix(coefficients(model, s = lambda)) # NOT doing [-1, ] to remove intercept!
    # NOTE: Even if there is a single model fitted, we can still request the coefficients() method validly with the s = lambda argument!
  }
  if(!intercept) {
    beta <- as.matrix(coefficients(model, s = lambda))[-1, ] # DOING [-1, ], note that coefficients will always provide an intercept...
  }
  
  ## Extracting the active set
  active_set <- which(beta != 0)
  
  # Taking intercept into account
  if(intercept) {
    active_set_noint <- active_set[-1]
  }
  if(!intercept) {
    active_set_noint <- active_set  # Then it's just the active set itself
  }
  
  # If the active set if not empty:
  if(length(active_set_noint) != 0) {
    # Constructing the matrix with only active columns X_active
    X_active <- X[, active_set]
    X_active_noint <- X[, active_set_noint]
    
    ## Producing predictions
    # Making specific matrix for predictions
    if(intercept) {
      X_noint <- as.matrix(X[, -1]) # Not necessary for predict() function
    }
    if(!intercept) {
      X_noint <- as.matrix(X)
    }
    
    # Predicting
    p_hat <- as.numeric(predict(object = model,
                                newx = X_noint,
                                type = "response",
                                s = lambda)) # NEW # TODO REMOVE (TEMP)
    # Forcing as.numeric() because diag(p_hat * (1 - p_hat)) not as nice with matrices
    
    ### Calculating degrees of freedom
    ## If alpha = 1 (LASSO): Trivial
    if(alpha == 1){
      df <- sum(beta != 0) # Amount of nonzero coefficients
    } 
    ## If alpha != 0: Calculate
    if(alpha != 1){
      ## Calculating Hat matrix (H)
      # Weight matrix W (diagonal), see GLM theory!
      W <- diag(p_hat * (1 - p_hat))
      
      # Calculating Hat (H) Matrix
      # NEW: Trying to fix wrong regularization
      I <- diag(length(active_set_noint))
      #H <- sqrt(W) %*% X_active %*% solve((t(X_active) %*% W %*% X_active)  + (lambda * (1 - alpha)/2) * diag(length(active_set))) %*% t(X_active) %*% sqrt(W)
      H <- sqrt(W) %*% X_active_noint %*% solve((t(X_active_noint) %*% W %*% X_active_noint)  + (lambda * (1 - alpha)/2) * I) %*% t(X_active_noint) %*% sqrt(W)  # NEW
      
      
      # Formula: W^(1/2) X_a (X'_a W X_a + lambda_2 * I)^-1 X'_a W^(1/2) (Yes transpose W, but it's diagonal, so it doesn't matter)
      # Note: lambda_2 = lambda (1 - alpha)/2 here
      
      # Degrees of freedom = trace of H
      H_trace <- sum(diag(H))
      
      # ADDING THE INTERCEPT AGAIN (UNREGULARIZED SO 1 EXTRA DF)
      if(intercept) {
        df <- 1 + H_trace
      }
      if(!intercept) {
        df <- H_trace
      }
    }
  }
  
  # If the active set length (without intercept) == 0
  if(length(active_set_noint) == 0) {
    if(intercept) {
      df <- 1
    } else if(!intercept){
      df <- 0
    }
  }
  
  # OUTPUT
  return(df)
}



