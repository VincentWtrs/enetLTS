lastbestindex_function <- function(zz, x, y, family, h, hsize, alpha, lambda, H2) {
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
  
  # while-loop keep running as long as objective function is bigger than tolerance unless max amount of csteps is reached
  countloop <- 0
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
    # Remember that a C-step also entails fitting a glmnet model!
  }
  return(list(lastindex = newindex, 
              objbest = objbest, 
              countloop = countloop, 
              residu = cstep.mod$residu, 
              beta = beta))
  
}