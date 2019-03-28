InitialSubset <- function(x, y, family, h, hsize, alpha, lambda, nsamp, scal, para, ncores, seed){
   # gives initial 500 subsamples after Two C Steps
   if (!is.null(seed)) set.seed(seed)
   if (family=="binomial"){
      index.subsets <- replicate(nsamp,c(sample(which(y==1),2),sample(which(y==0),2)))
   } else if (family=="gaussian"){
      index.subsets <- replicate(nsamp,sample.int(nrow(x), 3))
   }
   twoCstep <- function(c,x,y,family,h,hsize,alpha,lambda){
      ## C step 1
      if (family=="binomial"){
          Cstep1 <- CStep(x,y,family,index.subsets[,c],h,hsize,alpha,lambda/4,scal=FALSE)
      } else if (family=="gaussian"){
          Cstep1 <- CStep(x,y,family,index.subsets[,c],h,hsize,alpha,lambda/3,scal=FALSE)
      }
      indx1 <- Cstep1$index
      object1 <- Cstep1$object
      ## C step 2
      Cstep2 <- CStep(x,y,family,indx1,h,hsize,alpha,lambda/h,scal) # h observations
      indx2 <- Cstep2$index
      object <- Cstep2$object
      return(list(obj=object,indx=indx2))
   }
   if (para){
      subsets <- mclapply(1:nsamp,
                          FUN = twoCstep,
                          x = x, y = y,
                          family = family,
                          h = h, hsize = hsize,
                          alpha = alpha,
                          lambda = lambda,
                          mc.cores = ncores)
   } else {
      subsets <- lapply(1:nsamp,
                        FUN = twoCstep,
                        x = x, y = y,
                        family = family,
                        h = h, hsize = hsize,
                        alpha = alpha,
                        lambda = lambda)
   }
   return(list(subsets=subsets,index.subsets=index.subsets))
}

