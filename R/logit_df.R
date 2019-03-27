logit_df <- function(model, X, alpha = NULL, intercept){
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
  
  # Trying to retrieve alpha if missing in function call
  if(missing(alpha)){
    if(is.numeric(model$call$alpha)){
      alpha <- model$call$alpha # Trying to extract from the function call if possible
    } else { # But if it's another variable, it will return e.g. "myalpha" and will not be of any use
      stop("No alpha value given for logit_df(), also no value could be succesfully extracted from the model object, supply a valid alpha: [0, 1]")
    }
  }
  
  # TO DO: CHANGE BY CODE BELOW (model.matrix())
  # If the model was estimated using an intercept, the matrix calculations need to have a 1 in the column as well!
  if(intercept == TRUE){
    X <- cbind(1, X)
    colnames(X)[1] <- "X0" # Needs to be colnames because matrix!
  }
  
  # ALTERNATIVE TO INTERCEPT HASSLE (TO DO: TEST)
  #X <- model.matrix(model)

  # Extracting fitted beta (coefficients)
  if(intercept == TRUE){
    coefs <- as.matrix(coefficients(model))
  }
  if(intercept == FALSE){
    coefs <- as.matrix(coefficients(model))[-1, ] # EXCLUDING intercept, as.matrix() to get rid of sparse matrix formatting
    # Alternatively to get rid of intercept: model$beta
  }

  # Extracting the active subset
  active_set <- which(coefs != 0) # Nonzero coefs INDICES (Note they start from 1, not 0)
  X_active <- X[, active_set] # X subsetted with above indices
  
  # Extracting predictions
  p_hat <- as.numeric(predict(object = model,
                              newx = X,
                              type = "response")) 
  # Forcing as.numeric() because diag(p_hat * (1 - p_hat)) not as nice with matrices
  
  ### Calculating degrees of freedom
  ## If alpha = 1 (LASSO): Trivial
  if(alpha == 1){
    df <- model$df # TO DO: check with intercept, make sure it's consistent
  } 
  ## If alpha != 0: Calculate
  if(alpha != 1){
  ## Calculating Hat matrix (H)
  
  # Getting lambda
  lambda <- model$lambda
    
  # Weight matrix W (diagonal), see GLM theory!
  W <- diag(p_hat * (1 - p_hat))
  
  # Calculating Hat (H) Matrix
  H <- sqrt(W) %*% X_active %*% solve((t(X_active) %*% W %*% X_active)  + (lambda * (1 - alpha)/2) * diag(length(active_set))) %*% t(X_active) %*% sqrt(W)
  # Formula: W^(1/2) X_a (X'_a W X_a + lambda_2 * I)^-1 X'_a W^(1/2) (Yes transpose W, but it's diagonal, so it doesn't matter)
  # Note: lambda_2 = lambda (1 - alpha)/2 here
  
  # Degrees of freedom = trace of H
  H_trace <- sum(diag(H))
  df <- H_trace
  }
  
  # OUTPUT
  return(df)
}
   

  

# ##################################### OLD AND QUITE POSSIBLY WRONG BECAUSE ONE WORKS WITH GLMNET IN THE FITTING PROCEDURE NOT WITH ENETLTS ########
#   ### Extracting some parameters
#   ## enetLTS models
#   if(class(model) == "enetLTS"){
#     # Predictions
#     p_hat <- model$fitted.values # These are the reweigted fitted values!
#     
#     # Coefficients
#     coefs <- coefficients(model) # Gets coefficients (INCLUDING INTERCEPT)
#     coefs_nonzero <- coefs[coefs != 0]
#     index_coefs_nonzero <- as.numeric(names(coefs_nonzero)) # COLUMN NUMBERS of nonzero coefs
#     # Note the column numbers start from 1! and e.g. go to 51
#     
#     ## Getting subset out
#     wt <- model$raw.wt == 1
#     
#     X <- X[wt, ]
#     p_hat <- p_hat[wt]
#     
#     
#     
#     # Alpha, Lambda(reweighted)
#     lambda <- model$lambdaw
#     alpha <- model$alpha
#   } else {  ## Other kind of models
#     lambda <- model$lambda # UNDER THE ASSUMPTION ONLY A SINGLE LAMBDA IS GIVEN
#     p_hat <- predict(object = model, type = "response", newx = X[, -1]) # Without intercept fitting
#     
#   }
#   
#   # Extracting the active set
#   X <- X[, index_coefs_nonzero] # WITH OR WITHOUT INTERCEPT?!
#   
#   # Calculating W matrix
#   W <- diag(p_hat * (1 - p_hat))
#   
#   # Amount of obs
#   nobs <- length(p_hat)
#   
#   # Calculating Hat (H) Matrix
#   H <- sqrt(W) %*% X %*% solve((t(X) %*% W %*% X)  + lambda * (1 - alpha)/2) %*% t(X) %*% sqrt(W)
#   
#   # Calculating df 
#   H_trace <- sum(diag(H))
  

