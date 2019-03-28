uptrimMSE<- function(x, trim = 0.1){
  # computes trim% upper trimmed mean
  return(mean(x[x<quantile(x,1-trim)]))
}