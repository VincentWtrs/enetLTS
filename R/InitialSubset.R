InitialSubset <- function(x, y, family, h, hsize, alpha, lambda, nsamp, scal, para, ncores, seed){
   
   # Gives initial 500 subsamples after Two C Steps
   
   # Seed
   if (!is.null(seed)) {
      set.seed(seed)
   }
   
   ## Sampling observations
   # Case: Binomial sample 2 for each outcome
   if (family == "binomial") {
      index.subsets <- replicate(nsamp, c(sample(which(y == 1), 2), sample(which(y == 0), 2))) # Repeating it nsamp times
      
   # Case: Gaussian
   } else if (family == "gaussian") {
      index.subsets <- replicate(nsamp, sample.int(nrow(x), 3))  # Repeating it nsamp times
   }
   
   ## NEW: BEFORE THE twoCstep() FUNCTION WAS DEFINED HERE
   
   # Case: Parallel
   if (para) {
      subsets <- mclapply(1:nsamp, # Looping for nsamp times (Default: 500)
                          FUN = twoCstep, # This is the function to repeat
                          x = x, # From this point on some fixed parameters are provided
                          y = y,
                          family = family,
                          h = h, 
                          hsize = hsize,
                          alpha = alpha, # This will be the first alpha / lambda
                          lambda = lambda,
                          mc.cores = ncores,
                          index.subsets = index.subsets) # NEW: Added this because since we separated the function it now requires this as an input!
   # Case: Non-parallel
   } else {
      subsets <- lapply(1:nsamp, # Looping for nsamp times (Default: 500)
                        FUN = twoCstep, #  # This is the function to repeat
                        x = x, 
                        y = y,
                        family = family,
                        h = h,
                        hsize = hsize,
                        alpha = alpha,
                        lambda = lambda,
                        index.subsets = index.subsets) # NEW: Added this because since we separated the function it now requires this as an input!
   }
   
   # OUTPUT
   return(list(subsets = subsets,
               index.subsets = index.subsets))
}

