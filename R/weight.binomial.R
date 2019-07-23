weight.binomial <- function(x, y, beta, intercept, del){
  # TODO REMOVE DEBUG
  print(dim(cbind(1, x)))
  print(beta)
  print(head(x))
  
  
  if(intercept) {
    pi <- exp(cbind(1, x) %*% beta)/(1 + exp(cbind(1, x) %*% beta))
    res <- (y - pi) / sqrt(pi*(1-pi))
  } else{
    pi <- exp(x %*% beta[-1])/(1 + exp(x %*% beta[-1]))
    res <- (y - pi) / sqrt(pi * (1 - pi))
  }
  we <- as.integer(abs(res) <= qnorm(1 - del))
  return(we)
}