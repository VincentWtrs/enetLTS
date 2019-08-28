# COPY OF logit_df() TO EXPERIMENT WITH (COPY FROM 28/08/2019)
logit_df <- function(model, X, alpha = NULL, lambda = NULL, lasso_shortcut = TRUE, brute_force = FALSE, verbose = TRUE){
  ### logit_df() FUNCTION: calculating the effective degrees of freedom in a potentially penalized logistic regression using hat matrix trace
  
  ## INPUTS
  # model: A glmnet() model
  # X: the data matrix WITHOUT INTERCEPT
  # alpha: Alpha value used to fit the model (single value)
  # intercept: If intercept was used to fit the logit-glmnet model: set to TRUE
  # lambda: A specific lambda parameter to be given in case there is a model path supplied
  # lasso_shortcut: Uses unbiased estimator for the df (see Zhou & Hastie Talk 2005)
  # brute_force: Brute force: uses full expression with matrix inevrsion to get full H matrix. 
  
  ## OUTPUT
  # df: Estimated effective degrees of freedom
  
  ## NOTES
  # 1. lamba value can be extracted from glmnet model, alpha cannot in most circumstances
  # 2. The case for including an intercept is not so clean as it seems as it would require different
  # 3.   inner workings for the predict() function as it does not work with matrices with intercepts... might do this later
  # 3. I was afraid this might break down at the matrix inversion when |Active set| > n, but seems to be okay!
  
  ## TODO
  # 1. Extend the function to glm_df to open up options for Poisson e.g. 
  
  ## Defining detect intercept function
  detect_intercept <- function(model) {
    if(any(model$a0 !=  0)) {
      intercept <- TRUE
    } else {
      intercept <- FALSE
    }
    # OUTPUT
    return(intercept)
  }
  
  ## Defining detect regularization path function
  detect_regularization_path <- function(model) {
    if (length(model$lambda) > 1) {
      path <- TRUE
    } else {
      path <- FALSE
    }
  }
  
  ## Error handling
  # Checking if the model was fitted with an intercept
  intercept <- detect_intercept(model)
  
  # Checking if the model was fitted over a regularization path
  path <- detect_regularization_path(model)
  
  # Checking regularization path consistency with supplied lambda
  if (isTRUE(path)) {
    if (is.null(lambda)) {
      stop("logit_df was supplied a glmnet model fitted over a regularization path but no specific lambda was given to calculate the degrees of freedom")
    } else {
      # Setting the supplied lambda as the lambda to further use
      if(length(s) == 1){ # Additional check, we can only have a single lambda here
        lambda <- s
      }
    }
  }
  # TODO MAYBE THIS SHOULD BE MOVED IN WITH THE DETECT INTERCEPT FUNCTION
  
  
  # Trying to retrieve alpha if missing in function call
  if (missing(alpha)) {
    if(is.numeric(model$call$alpha)){ # If the function was called with explicit numeric value it can retrieve a meaningfull alpha
      alpha <- model$call$alpha # Trying to extract from the function call if possible
    } else { # But if it's another variable, it will return e.g. "my_alpha" and will not be of any use
      stop("No alpha value given for logit_df(), also no value could be succesfully extracted from the model object, supply a valid alpha: [0, 1]")
    }
  }
  
  
  ## Construction of the design matrix
  # If there is an intercept, add a first column of 1
  if(isTRUE(intercept)) {
    X <- cbind(1, X)
    colnames(X)[1] <- "X0"
  } # If there is no intercept, the X matrix is good as is
  
  ## Extracting fitted beta coefficients
  # If intercept, get the beta0 in there as well!
  if (isTRUE(intercept)) {
    beta <- as.matrix(coefficients(model, s = lambda)) # NOT doing [-1, ] to remove intercept!
    # NOTE: Even if there is a single model fitted, we can still request the coefficients() method validly with the s = lambda argument!
  }
  if (isFALSE(intercept)) {
    beta <- as.matrix(coefficients(model, s = lambda))[-1, ] # DOING [-1, ], note that coefficients will always provide an intercept...
  }
  
  ## Extracting the active set
  active_set <- which(beta != 0)
  
  # Taking intercept into account
  if(isTRUE(intercept)) {
    active_set_noint <- active_set[-1] # Removing first element (intercept)
  }
  if(isFALSE(intercept)) {
    active_set_noint <- active_set  # Then it's just the active set itself
  }
  
  # If the active set if not empty:
  if(length(active_set_noint) != 0) {
    
    # Constructing the matrix with only active columns X_active
    X_active <- X[, active_set]
    X_active_noint <- X[, active_set_noint]
    
    ## Producing predictions
    # NOTE: Since predict() will handle the active set itself properly we just need to make sure we give the full matrix but WITHOUT CONSTANT COLUMN!
    if (isTRUE(intercept)) {
      X_noint <- as.matrix(X[, -1]) # Removing the constant column from the X matrix which is (1, X)
    }
    if (isFALSE(intercept)) {
      X_noint <- as.matrix(X)
    }
    
    # Predicting
    p_hat <- as.numeric(predict(object = model,
                                newx = X_noint, # It wil automatically add intercept if it was ever fitted with one!
                                type = "response",
                                s = lambda)) # Fitting a single lambda
    
    ### Calculating degrees of freedom
    ## If alpha = 1 (LASSO): Shortcut calculation
    if(isTRUE(lasso_shortcut)) {
      if(alpha == 1) { # Amount of nonzero coefficients
        df <- sum(beta != 0)
      }
    }
    
    ## If alpha != 0: Calculate using weighted hat matrix
    if (alpha != 1) {
      if (isTRUE(brute_force)) {
        ## Calculating Hat matrix (H)
        # Weight matrix W (diagonal), see GLM theory!
        W <- diag(p_hat * (1 - p_hat))
        W_sqrt <- diag(sqrt(p_hat * (1 - p_hat))) # Making square root beforehand (Potentially faster)
        I <- diag(length(active_set_noint))
        H <- W_sqrt %*% X_active_noint %*% solve((t(X_active_noint) %*% W %*% X_active_noint)  + (lambda * (1 - alpha)/2) * I) %*% t(X_active_noint) %*% W_sqrt  # NEW
        # NOTE: SINCE NO REGULARIZATION IS APPLIED ON THE INTERCEPT: EXCLUDE FROM DF CALCULATION (EXCEPT PREDS) AND ADD 1 AT THE END
        # NOTE: ONLY FOR GLM PHI = CONSTANT MODEL, OTHERWISE: ADD 1 EXTRA AS WELL
        # NOTE: lambda_2 = lambda (1 - alpha)/2 here
        
        H_trace <- sum(diag(H))
        
      } else {
        print('Using Cholesky decomposition!')
        
        # Defining intermediate matrices
        w_ii <- sqrt(p_hat * (1 - p_hat)) # Square root Weight vector
        W <- diag(p_hat * (1 - p_hat)) # Weight matrix
        X_active_noint_tilde <- w_ii * X_active_noint # X_tilde in most literature
        I <- diag(length(active_set_noint)) # Identity matrix of size of active set
        
        # Choleski Decomposition
        A <- t(X_active_noint) %*% W %*% X_active_noint + (lambda * (1 - alpha)/2)  # Matrix that needs would need to be inverted 'A'
        U <- base::chol.default(A)
        Z <- base::forwardsolve(t(U), t(X_active_noint_tilde))
        
        # Column Sums to get diagonal matrix
        diag_ZZ <- colSums(Z^2)
        
        # To become trace: Sum
        H_trace <- sum(diag_ZZ) # It's the same as H_trace now!
      }
      
      # ADDING THE INTERCEPT AGAIN (UNREGULARIZED SO 1 EXTRA DF)
      if(isTRUE(intercept)) {
        df <- 1 + H_trace
      }
      if(isFALSE(intercept)) {
        df <- H_trace
      }
    }
  }
  
  # If the active set length (without intercept) == 0
  if(length(active_set_noint) == 0) {
    if (intercept) {
      df <- 1
    } else if(!intercept){
      df <- 0
    }
  }
  
  # OUTPUT
  print(df)
  print(paste0("Nonzeros:", model$df))
  return(df)
}



