winsorize.default <- function(x, standardized = FALSE, centerFun = median,
                              scaleFun = mad, const = 2,
                              return = c("data", "weights"), ...){
  ## initializations
  standardized <- isTRUE(standardized)
  if(standardized) return <- match.arg(return)
  else {
    # standardize data
    x <- robStandardize(x, centerFun=centerFun, scaleFun=scaleFun, ...)
    center <- attr(x, "center")
    scale <- attr(x, "scale")
  }
  ## winsorize standardized data
  #   ind <- abs(x) > const           # observations in 'x' that need to be shrunken
  #   x[ind] <- const * sign(x[ind])  # winsorize
  weights <- pmin(const / abs(x), 1)
  if(standardized && return == "weights") return(weights)
  x <- weights * x
  ## finalizations
  if(!standardized) {
    # transform back to original scale and remove attributes
    x <- c(x * scale + center)
  }
  x
}