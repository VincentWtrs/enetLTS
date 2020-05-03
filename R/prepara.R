
prepara <- function(x, y, family, index = NULL, robu = NULL){
  
  ## GOAL
  # Handles centering and scaling of the input data using classic or robust estimators of location and scale.
  # ... Default: Classic scaling. If robu = 1: robust scaling is used based on median for y (Gaussian only)
  # ... and median + MAD for the predictors.
  
  ## INPUTS
  # x: Predictor matrix X
  # y: Outcome matrix y
  # family: GLM family argument, currently supported: "gaussian", "binomial"
  # index: Indices to apply on the standardization of the data
  # robu: Use of robust scale/location standardization (1) or not (0)
  
  ## STANDING ISSUES:
  # THIS THING FAILS ON DUMMY VARS, BECAUSE THE MAD MIGHT BE EQUAL TO 0. I THINK THAT -FOR THE DUMMY COLUMNS- WE SHOULD JUST USE ORDINARY SCALING
  # ADDITIONAL MESS IS: KEEPING THE STUFF IN THE RIGHT ORDER (COLUMN-WISE...)
  
  original_col_order <- colnames(x)  # Necessary to put all result in their original order!
  
  factor_handling <- function(x) {
    # Count unique values per column
    col_nunique <- apply(X=x, MARGIN=2, function(x) length(unique(x)))
    
    # Checking if there are any potential factor columns by checking amount of cols that have 2 values only
    amount_factors <- length(col_nunique[col_nunique==2])
    
    # If there are any: they need special treatment:
    if (amount_factors > 0) {
      # Gathering the column positions of those factors
      factor_col_indices <- which(col_nunique == 2)
      
      return(factor_col_indices)  # Return these (No return if no factors)
    }
  }
  
  # Checking if there are factors
  factor_col_indices <- factor_handling(x)  # If there are no factor cols, return will be NULL
  if (!is.null(factor_col_indices)) {
    x_factor <- as.matrix(x[, factor_col_indices])  # Forcing as.matrix() s.t. apply() later on keeps working
    x <- as.matrix(x[, -factor_col_indices])  # Forcing as.matrix() s.t. apply() later on keeps working
  }

  # If no robu parameter given: will go to default: NULL: set robu <- 0
  if (is.null(robu)) {
    robu <- 0
  }
  
  # 
  if (is.null(index)) {
    
    if (robu > 0){
      if (family == "binomial") {
        muy <- y
      } else if (family == "gaussian") {
        muy <- median(y)
      }
      mux <- apply(x, 2, median)   ## to go back original coef
      sigx <- apply(x, 2, mad) 
      if (!is.null(factor_col_indices)) {
        mux_factor <- apply(x_factor, 2, mean)
        mux <- c(mux, mux_factor)
        mux <- mux[original_col_order]  # Can be done with named vectors
        
        sigx_factor <- apply(x_factor, 2, sd)
        sigx <- c(sigx, sigx_factor)
        sigx <- sigx[original_col_order]  # Can be done with named vectors
        
        x <- cbind(x, x_factor)
        x <- x[, original_col_order]

      }
      
    } else {
      if(family == "binomial"){
        muy <- y
      } else if(family =="gaussian") {
        muy <- mean(y)
      }
      mux <- apply(x, 2, mean)   ## to go back original coef
      sigx <- apply(x, 2, sd)
    }
  } else { # If index not NULL:
    if (robu > 0){
      if (family == "binomial") {
        muy <- y
      } else if (family == "gaussian") {
        muy <- median(y[index])
      }
      mux <- apply(x[index, ], 2, median)
      sigx <- apply(x[index, ], 2, mad)
      if (!is.null(factor_col_indices)) {
        mux_factor <- apply(x_factor[index, ], 2, mean)
        mux <- c(mux, mux_factor)
        mux <- mux[original_col_order]
        
        sigx_factor <- apply(x_factor[index, ], 2, sd)
        sigx <- c(sigx, sigx_factor)
        sigx <- sigx[original_col_order]
        
        x <- cbind(x, x_factor)
        x <- x[, original_col_order]
      }
    } else {
      if (family == "binomial") {
        muy <- y
      } else if (family == "gaussian") {
        muy <- mean(y[index])
      }
      mux <- apply(x[index, ], 2, mean)
      sigx <- apply(x[index, ], 2, sd)
    }
  }
  # Effectively scaling X
  xnor <- scale(x, mux, sigx)
  # NOTE: x is being scaled as a whole, so the returning X has the same sample size as the input X!
  
  if (family == "binomial") {
    ycen <- y
    } else if (family == "gaussian") {
    ycen <- scale(y, muy, FALSE)
  }
  
  # OUTPUT
  return(list(xnor = xnor,
              ycen = ycen, 
              mux = mux, 
              sigx = sigx,
              muy = muy))
}