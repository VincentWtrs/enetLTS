prepara <- function(x, y, family, index = NULL, robu = NULL){
  # default clasical scale, robust scale for robu=1

  if (is.null(robu)) robu=0
  if (is.null(index)) {
    if (robu>0){
      if(family=="binomial"){
        muy=y
      } else if(family=="gaussian"){
        muy <- median(y)
      }
      mux <- apply(x,2,median)   ## to go back original coef
      sigx <- apply(x,2,mad)
    } else {
      if(family=="binomial"){
        muy=y
      } else if(family=="gaussian"){
        muy <- mean(y)
      }
      mux <- apply(x,2,mean)   ## to go back original coef
      sigx <- apply(x,2,sd)
    }
  } else {
    if (robu>0){
      if(family=="binomial"){
        muy=y
      } else if(family=="gaussian"){
        muy <- median(y[index])
      }
      mux <- apply(x[index,],2,median)
      sigx <- apply(x[index,],2,mad)
    } else {
      if(family=="binomial"){
        muy=y
      } else if(family=="gaussian"){
        muy <- mean(y[index])
      }
      mux <- apply(x[index,],2,mean)
      sigx <- apply(x[index,],2,sd)
    }
  }
  xnor <- scale(x,mux,sigx)
  if(family=="binomial"){
    ycen <- y}
  else if(family=="gaussian"){
    ycen <- scale(y,muy,FALSE)}
  return(list(xnor=xnor,ycen=ycen,mux=mux,sigx=sigx,muy=muy))
}