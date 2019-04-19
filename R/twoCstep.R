twoCstep <- function(c, x, y, family, h, hsize, alpha, lambda, index.subsets, scal) {
  # NEW: needed to add the index.subsets as an input because since separation the function could otherwise not see this input
  
  ## C-step 1
  # Case: Binomial
  if (family=="binomial") {
    Cstep1 <- CStep(x = x,
                    y = y,
                    family = family,
                    indx = index.subsets[, c],
                    h = h,
                    hsize = hsize,
                    alpha = alpha,
                    lambda = lambda/4, # 4 (2 + 2) observations Binomial case
                    scal = FALSE)
    
  # Case: Gaussian
  } else if (family=="gaussian") {
    Cstep1 <- CStep(x = x,
                    y = y,
                    family = family,
                    indx = index.subsets[, c],
                    h = h,
                    hsize = hsize,
                    alpha = alpha,
                    lambda = lambda/3, # 3 observations enough in Gaussian case
                    scal = FALSE)
  }
  
  # Gathering results from C-step 1
  indx1 <- Cstep1$index
  object1 <- Cstep1$object
  
  ## C-step 2 (Same for both cases)
  Cstep2 <- CStep(x = x,
                  y = y,
                  family = family,
                  indx = indx1,
                  h = h,
                  hsize = hsize,
                  alpha = alpha,
                  lambda = lambda/h, # h observations
                  scal = scal) 
  
  # Gathering results from C-step 2
  indx2 <- Cstep2$index
  object <- Cstep2$object
  
  return(list(obj = object, indx = indx2))
}