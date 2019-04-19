selectbest10 <- function(x, y, family, h, hsize, alpha, lambda, nsamp, s1, para, ncores, scal, seed){
  obj <- NULL
  all_subsets <- InitialSubset(x,y,family,h,hsize,alpha,lambda,nsamp,para,ncores,scal,seed)
  subsets <- all_subsets$subsets
  index.subsets <- all_subsets$index.subsets
  if(para){
    obj <- unlist(mclapply(1:nsamp, function(ob,sub){
      ob_val <- subsets[[ob]]$obj
    }, subsets, mc.cores = ncores, mc.allow.recursive = FALSE))
  }else{
    for (i in 1:nsamp){ obj <- c(obj,subsets[[i]]$obj) }
  }
  if(family == "binomial"){
    obj_sorted <- sort(obj,decreasing=TRUE,index.return=TRUE)
  } else if(family=="gaussian"){
    obj_sorted <- sort(obj,decreasing=FALSE,index.return=TRUE)
  }
  obj <- obj_sorted$x[1:s1]
  s1_new <- length(obj[!is.infinite(obj)])
  idx <- obj_sorted$ix[1:s1_new]
  if(s1_new == 0){
    stop(paste("Model is not suitable for alpha",alpha,"lambda",lambda,"for this data set. Choose another lambda."))
  }
  if(para){
    bestindex <- mclapply(1:s1_new, function(c,idx,subsets) {
      indx <- subsets[[idx[c]]]$indx
    },idx,subsets,mc.cores = ncores)
  }else{
    bestindex <- lapply(1:s1_new, function(c,idx,subsets) {
      indx <- subsets[[idx[c]]]$indx
    },idx,subsets)
  }
  return(list(idxbest = bestindex, 
              s1_new = s1_new,
              subsets = subsets,
              index.subsets = index.subsets))
}