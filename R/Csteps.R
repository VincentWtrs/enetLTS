beginningCstep_OLD <- function(x, y, family, h, hsize, alpha, lambda, nsamp, s1, ncores, csteps, tol, scal, para, seed){
 
 ## internal function for Cstep and warmCsteps
 #  source("objectiveFunc.R")
 #  source("InitialSubsets.R")
 #  source("Csteps.R")
 #  source("utilities.R")
   
   H2 <- selectbest10(x,y,family,h,hsize,alpha,lambda,nsamp,s1,para,ncores,scal,seed) 
   if (para){ 
      lastbestindex <- mclapply(1:s1, function(zz,x,y,family,h,hsize,alpha,lambda,H2) {
         indexsubbest <- H2$idxbest[[zz]]
         objbest <- tol
         cstep.mod <- CStep(x,y,family,indexsubbest,h,hsize,alpha,lambda/h,scal)
         countloop <- 0
         while ((cstep.mod$object>objbest) & (countloop<csteps)){ 
            countloop <- countloop+1
            objbest <- cstep.mod$object 
            newindex <- cstep.mod$index  
            beta <- cstep.mod$beta
            cstep.mod <- CStep(x,y,family,newindex,h,hsize,alpha,lambda/h,scal)
         }
         return(list(lastindex=newindex,objbest=objbest,countloop=countloop,
                     residu=cstep.mod$residu,beta=beta))
      },x,y,family,h,hsize,alpha,lambda,H2,mc.cores = ncores) 
   }else{ # not parallel
      lastbestindex <- lapply(1:s1, function(zz,x,y,family,h,hsize,alpha,lambda,H2) {
         indexsubbest <- H2$idxbest[[zz]] 
         objbest <- tol
         cstep.mod <- CStep(x,y,family,indexsubbest,h,hsize,alpha,lambda/h,scal)
         countloop <- 0
         while ((cstep.mod$object>objbest) & (countloop<csteps)){
            countloop <- countloop+1
            objbest <- cstep.mod$object 
            newindex <- cstep.mod$index  
            beta <- cstep.mod$beta
            cstep.mod <- CStep(x,y,family,newindex,h,hsize,alpha,lambda/h,scal)
         }
         return(list(lastindex=newindex,objbest=objbest,countloop=countloop,
                     residu=cstep.mod$residu,beta=beta))
      },x,y,family,h,hsize,alpha,lambda,H2)
   } 
   obj <- NULL
   for (i in 1:s1){
      obj <- c(obj,lastbestindex[[i]]$objbest)
   }
   whichbestindex <- sort(obj,decreasing=TRUE,index.return=TRUE)$ix[1]
   index <- lastbestindex[[whichbestindex]]$lastindex
   resid <- lastbestindex[[whichbestindex]]$residu
   # beta <- lastbestindex[[whichbestindex]]$beta
   return(list(index=index,resid=drop(resid))) #,s1=s1,beta=beta))
}