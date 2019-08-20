weight.binomial <- function(x, y, beta, intercept, del) {
  
  ## This function gives 0/1 weights based on the Pearson residuals, which are assumed to be approximately standard normal.S
  
  ## Case: Intercept (Always for logistic regression)
  if(intercept) {
    # Calculating predicted probabilities
    pi <- exp(cbind(1, x) %*% beta) / (1 + exp(cbind(1, x) %*% beta))
    
    # Calculating Pearson Residual (Non-studentized)
    res <- (y - pi) / sqrt(pi*(1-pi))
    
    ## Catching some NAs for very high predicted probabilities
    if (sum(is.na(res)) > 1){
      # Catching the indices that are NA
      res_na_indx <- which(is.na(res))
      
      # Catching indices that are very sure of their predictions and are CORRECT
      correct_sure_indx <- which((y - pi) < 0.000000001)
      
      # Intersection of these both
      intersection_na_sure <- intersect(res_na_indx, correct_sure_indx)
      
      # Assigning 0 to these (we are very sure and correct!)
      res[intersection_na_sure] <- 0
      
      # Taking difference with the sure ones
      incorrect_sure_indx <- setdiff(res_na_indx, correct_sure_indx)
      
      # Setting those incorrect sure ones rather big
      print("The residuals before replacing had the following y - pi_hat values:")
      print("printing y[incorrect_sure_indx]")
      print(y[incorrect_sure_indx])
      print("Printing pi[incorrect_sure_indx]")
      print(pi[incorrect_sure_indx])
      
      print(y[incorrect_sure_indx] - pi[incorrect_sure_indx])
      res[incorrect_sure_indx] <- 3 # 3 Standard deviations away so will be flagged
    }

  ## Case: No intercept 
  } else {
    
    # Calculating predicted probabilities
    pi <- exp(x %*% beta[-1])/(1 + exp(x %*% beta[-1]))
    res <- (y - pi) / sqrt(pi * (1 - pi))
  }
  
  
  we <- as.integer(abs(res) <= qnorm(1 - del))
  
  return(we)
}