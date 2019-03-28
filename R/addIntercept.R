## add intercept column to design matrix
addIntercept <- function(x, check = FALSE) {
  if(!check || is.na(match("(Intercept)", colnames(x)))) {
    cbind("(Intercept)"=rep.int(1, nrow(x)), x)
  } else x
}