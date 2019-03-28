nonzeroCoef.enetLTS  <- function(beta){
  beta <- as.matrix(beta)
  beta <- abs(beta)>0      # this is sparse
  beta <- which(beta)
  names(beta) <- 1:length(beta)
  return(beta)
}

