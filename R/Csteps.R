beginningCstep <- function(x, y, family, h, hsize, alpha, lambda, nsamp, s1, ncores, csteps, tol, scal, para, seed){
 
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

#####################################################################################################################

CStep <- 
function(x,y,family,indx,h,hsize,alpha,lambda,scal)
{
   ## internal function

   # require(glmnet)
   # source("utilities.R")
   # source("objectiveFunc.R")
   
   n <- nrow(x)
   if (scal){
      scl <- prepara(x,y,family,indx,robu=0)
      xs <- scl$xnor
      ys <- scl$ycen
      if (family=="binomial"){
         fit <- glmnet(xs[indx,],ys[indx],family,alpha=alpha,lambda=lambda,standardize=FALSE,intercept=FALSE)
         beta <- matrix(fit$beta)
         resid <- -(ys * xs %*% beta) + log(1+exp(xs %*% beta))
         if(all(beta==0)){return(list(object=-Inf,index=indx,residu=resid,beta=beta))} 
         resid.sort <- sort(resid,decreasing=FALSE,index.return=TRUE) 
         h0 <- floor((length(y[y==0])+1)*hsize)
         h1 <- h-h0
         index0 <- resid.sort$ix[y[resid.sort$ix]==0][1:h0]
         index1 <- resid.sort$ix[y[resid.sort$ix]==1][1:h1]
         indxnew <- c(index0,index1)
      }else if(family=="gaussian"){
         fit <- glmnet(xs[indx,],ys[indx],family,alpha=alpha,lambda=lambda,standardize=FALSE,intercept=FALSE)
         beta <- matrix(fit$beta)
         resid <- ys - predict(fit,xs,exact=TRUE)
         resid.sort <- sort(abs(resid),index.return=TRUE)
         indxnew <- resid.sort$ix[1:h]
      }
      obj <- Objval(xs,ys,family,beta,indxnew,alpha,lambda)
   }else{
      if (family=="binomial"){
         fit <- glmnet(x[indx,],y[indx],family,alpha=alpha,lambda=lambda,standardize=FALSE,intercept=FALSE)
         beta <- matrix(fit$beta)
         resid <- -(y * x %*% beta) + log(1+exp(x %*% beta))
         if(all(beta==0)){return(list(object=-Inf,index=indx,residu=resid,beta=beta))}
         resid.sort <- sort(resid,decreasing=FALSE,index.return=TRUE) 
         h0 <- floor((length(y[y==0])+1)*hsize)
         h1 <- h-h0
         index0 <- resid.sort$ix[y[resid.sort$ix]==0][1:h0]
         index1 <- resid.sort$ix[y[resid.sort$ix]==1][1:h1]
         indxnew <- c(index0,index1)
      }else if(family=="gaussian"){
         fit <- glmnet(x[indx,],y[indx],family,alpha=alpha,lambda=lambda,standardize=FALSE,intercept=FALSE)
         beta <- matrix(fit$beta)
         resid <- y - predict(fit,x,exact=TRUE)
         resid.sort <- sort(abs(resid),index.return=TRUE)
         indxnew <- resid.sort$ix[1:h]
      }
      obj <- Objval(x,y,family,beta,indxnew,alpha,lambda)
   }
   return(list(object=obj,index=indxnew,residu=resid,beta=beta))
}





