coef.enetLTS <- function(object, vers=c("reweighted","raw"), zeros = TRUE, ...){
   vers <- match.arg(vers)
   nbeta <- predict.enetLTS(object, 
                            newX = object$inputs$xx, 
                            vers = vers, 
                            type = "coefficients", ...)
   nbeta <- as.numeric(unlist(nbeta))
   if (isTRUE(zeros)) { # TO DO ?
      nbeta <- nbeta
      names(nbeta) <- 0:length(nbeta) # Changing to start from X0
   } else if (!isTRUE(zeros)) { # TO DO: check
      namesbeta <- which(nbeta != 0)
      nbeta <- nbeta[nbeta != 0]
      #names(nbeta) <- namesbeta # Check if replacing this helps
      names(nbeta) <- 0:length(nbeta)
   }
   return(nbeta)
}


