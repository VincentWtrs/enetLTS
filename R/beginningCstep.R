beginningCstep <- function(x, y, family, h, hsize, alpha, lambda, nsamp, s1, ncores, csteps, tol, scal, para, seed){
  
  ### beginningCstep() FUNCTION
  ## This function is called by warmCsteps() with following inputs: x = x, y = y, family = family, h = h, hsize = hsize,
  # alpha = alpha, lambda = lambda, nsamp = nsamp, s1 = s1, ncores = ncores, csteps = csteps, tol = tol, scal = scal, para = para, seed = seed
  
  ## From the KHF (2017) paper: this beginning Cstep is only done for the first combination of alpha, lambda
  
  # Selects best 10 subsets
  H2 <- selectbest10(x = x, 
                     y = y, 
                     family = family, 
                     h = h, 
                     hsize = hsize, 
                     alpha = alpha, 
                     lambda = lambda, 
                     nsamp = nsamp, 
                     s1 = s1, 
                     para = para, 
                     ncores = ncores, 
                     scal = scal, 
                     seed = seed) # Seems to work 
  # Note: warning (lognet() .... < 8, dangerous grounds) are thrown because the fit is made with only 4 observations
  # H2: H2$idxbest: list of 10 (s1_new) with 75 indices; H2$s1_new: numeric stating s1_new (might be different form s1 if there were infs in objective)
  # H2$subsets: list of 500, each of the 500 elements is a list of 2 containing 75 indices and their objective functions
  # H2$index.subsets: one single matrix of 4 rows (size elemental subsets) and nsamp columns
  # names(H2): "idxbest"       "s1_new"        "subsets"       "index.subsets"
  
  ## H2 EXAMPLES:
  # H2$idxbest[[1]]: Contains 75 observation indices of the first of s1 (10) subsets, ... H2$idxbest[[10]]: Contains h (75) observation indices of the last of s1 (10) subsets
  # H2$subsets[[1]]: Contains $obj: with the objective function value (Bianco-Yohai) for subset 1, $indx: vector with 75 observation indices, ... , up to H2$subsets[[500]].
  # H2$index.subsets: Matrix (4 x nsamp (4 x 500)) that contains the observation indices of the (initial) elemental subsets (TO DO ?)
  
  # IT SEEMS THAT ONLY $idxbest IS REALLY USED FURTHER ON!
  
  ## Running a lot of additional C-steps on the 10 Best subsets!
  # Parallel Case (see futher for non-parallel case)
  if (para) {
    lastbestindex <- mclapply(1:s1, function(zz, x, y, family, h, hsize, alpha, lambda, H2) { ## DEFINNG FUNCTION TO BE RAN WITHIN MCLAPPLY
      indexsubbest <- H2$idxbest[[zz]]
      objbest <- tol
      cstep.mod <- CStep(x = x, 
                         y = y, 
                         family = family, 
                         indx = indexsubbest, 
                         h = h, 
                         hsize = hsize, 
                         alpha = alpha, 
                         lambda = lambda/h, 
                         scal = scal) # Running a C-step
      countloop <- 0 # Initiating countloop at 0
      
      # while loop to run a set amount of csteps or stop early when tolerance level is reached
      while ((cstep.mod$object > objbest) & (countloop < csteps)) {
        countloop <- countloop + 1
        objbest <- cstep.mod$object
        newindex <- cstep.mod$index
        beta <- cstep.mod$beta
        cstep.mod <- CStep(x = x, 
                           y = y, 
                           family = family, 
                           indx = newindex, 
                           h = h, 
                           hsize = hsize, 
                           alpha = alpha, 
                           lambda = lambda/h, 
                           scal = scal)
      }
      # OUPUT
      return(list(lastindex = newindex, 
                  objbest = objbest, 
                  countloop = countloop, 
                  residu = cstep.mod$residu, 
                  beta = beta))
    }, x, y, family, h, hsize, alpha, lambda, H2, mc.cores = ncores) ## END OF DEFINING FUNCTION WITHIN MCLAPPLY
    
  # Non-parallel case
  } else {
    # Looping over all of the s1 (10) best subsets
    lastbestindex <- lapply(1:s1, function(zz, x, y, family, 
                                           h, hsize, alpha, lambda, H2) { # Beginning function definition
      indexsubbest <- H2$idxbest[[zz]] # zz is the number of the 10 subsets (1, ..., 10)
      objbest <- tol
      cstep.mod <- CStep(x = x, 
                         y = y, 
                         family = family, 
                         indx = indexsubbest, 
                         h = h, 
                         hsize = hsize, 
                         alpha = alpha, 
                         lambda = lambda/h, 
                         scal = scal)
      countloop <- 0
      # while-loop keep running as long as objective function is bigger than tolerance unless max amount of csteps is reached
      while ((cstep.mod$object > objbest) & (countloop < csteps)) {
        countloop <- countloop + 1 # Need to keep track of runs to stop on time
        objbest <- cstep.mod$object # cstep.mod will be overwritten each time but no worries, the C-steps are decreasing the objective each time anyway!
        newindex <- cstep.mod$index # Keeping the indices of course
        beta <- cstep.mod$beta # And the fitted betas
        cstep.mod <- CStep(x = x, 
                           y = y, 
                           family = family, 
                           indx = newindex, 
                           h = h, 
                           hsize = hsize, 
                           alpha = alpha, 
                           lambda = lambda/h, # Correct scaling of lambda as usual
                           scal = scal)
        # Remember that a C-step also entails fitting the model!
      }
      return(list(lastindex = newindex, 
                  objbest = objbest, 
                  countloop = countloop, 
                  residu = cstep.mod$residu, 
                  beta = beta))
    # Ending function definition
    } , x, y, family, h, hsize, alpha, lambda, H2) # Passing params
  } # END else
  
  
  ### Getting the results
  # Getting objective functions for each of 10 subsets after all the many Csteps
  obj <- NULL
  for (i in 1:s1) { # For each of the 10 best subsets
    obj <- c(obj, lastbestindex[[i]]$objbest) # Taking objbest (defined as output see above) out of the result of the lapply loop
  }
  
  ## GETTING THE ABSOLUTE WINNER OUT OF THE 10
  whichbestindex <- sort(obj, decreasing = TRUE, index.return = TRUE)$ix[1] # This gives a single number
  
  # Getting indices and residuals
  index <- lastbestindex[[whichbestindex]]$lastindex # Gives h indices of a clean dataset to be used
  resid <- lastbestindex[[whichbestindex]]$residu # Gives the n (!!!!!!) residuals (technically in a matrix, will make it to vector)
  
  return(list(index = index, 
              resid = drop(resid))) # drop(): drops the dimensionalities if possible, hence makes the matrix a vector
}