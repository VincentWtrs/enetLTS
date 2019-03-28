weight.gaussian <- function(resi, ind, del){
  if(is.logical(ind)){
    h <- length(which(ind==TRUE))
  }else{
    h <- length(ind)
  }
  n <- length(resi)
  mu <- mean(resi[ind])
  rc <- (resi - mu)
  qn <- qnorm((h+n)/ (2*n))                         # required quantile
  cdelta <- 1 / sqrt(1 - (2*n)/(h/qn) * dnorm(qn))
  s <- sqrt(mean(rc[ind]^2)) * cdelta
  we <- as.integer(abs(rc/s) <= qnorm(1-del))
  out <- list(we=we,mu=mu,s=s)
  
  return(out)
}