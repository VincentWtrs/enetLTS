enetLTS_reweight_results <- function(xx, yy, family, h, indexbest, alphabest, lambdabest, lambdaw = NULL, del, intercept, scal, type_lambdaw) {
  ## WRAPPER FUNCTION
  
  ### STEP: REWEIGHTING AND REFITTING
  ## Running function
  # Case for missing lambdaw (Calling function without lambdaw (commented out))
  if(is.null(lambdaw)){
    enetLTS_reweighting_refit <- enetLTS_reweighting_refitting(xx = xx,
                                                               yy = yy,
                                                               family = family, 
                                                               indexbest = indexbest, 
                                                               alphabest = alphabest, 
                                                               lambdabest = lambdabest, 
                                                               #lambdaw = lambdaw, # If missing don't call with lambdaw
                                                               del = del,
                                                               intercept = intercept,
                                                               scal = scal,
                                                               type_lambdaw = type_lambdaw)
    # Case for non-missing lambdaw
  } else {
    enetLTS_reweighting_refit <- enetLTS_reweighting_refitting(xx = xx,
                                                               yy = yy,
                                                               family = family, 
                                                               indexbest = indexbest, 
                                                               alphabest = alphabest, 
                                                               lambdabest = lambdabest, 
                                                               lambdaw = lambdaw, # otherwise call with
                                                               del = del,
                                                               intercept = intercept,
                                                               scal = scal,
                                                               type_lambdaw = type_lambdaw)
  }
  
  ## Extracting results back into the enetLTS function
  # Common for both families (Gaussian, Binomial)
  raw.wt <- enetLTS_reweighting_refit$raw.wt
  wgt <- enetLTS_reweighting_refit$wgt
  a00 <- enetLTS_reweighting_refit$a00
  a0 <- enetLTS_reweighting_refit$a0
  coefficients <-  enetLTS_reweighting_refit$coefficients
  raw.coefficients <- enetLTS_reweighting_refit$raw.coefficients
  alphabest <- enetLTS_reweighting_refit$alphabest
  alphaw_best <- enetLTS_reweighting_refit$alphaw_best
  lambdaw_best <- enetLTS_reweighting_refit$lambdaw_best
  raw.residuals <- enetLTS_reweighting_refit$raw.residuals
  raw.rmse <- enetLTS_reweighting_refit$raw.residuals
  reweighted.residuals <- enetLTS_reweighting_refit$reweighted.residuals
  num.nonzerocoef <- enetLTS_reweighting_refit$num.nonzerocoef
  reweighted_cv <- enetLTS_reweighting_refit$reweighted_cv
  
  # For Gaussian only:
  if(family == "gaussian") {
    raw.rmse <- enetLTS_reweighting_refit$raw.rmse
    reweighted.rmse <- enetLTS_reweighting_refit.reweighted.rmse
  }
  
  ## Plotting CV Plot 
  plot(reweighted_cv)
  print(paste0("The optimal reweighted lambda is: ", lambdaw_best)) # Printing chosen lambdaw
  print(paste0("The optimal alpha used in the reweighting step is: ", alphaw_best)) # Printing chosen alpha
  # Note: I put this here at the end because otherwise I would need to repeat the statement multiple times
  
  ## PREPARING OUTPUT (Already given in the reweighting step...)
  # Counting number of nonzero coefficients
  #num.nonzerocoef <- sum(coefficients != 0)
  
  ## NEW FUNCTION TO HANDLE RESULTS!
  # Looping over iCS
  results <- enetLTS_results_handling(xx = xx,
                                      yy = yy,
                                      indexbest = indexbest,
                                      intercept = intercept,
                                      family = family,
                                      h = h,
                                      a0 = a0,
                                      a00 = a00,
                                      coefficients = coefficients,
                                      raw.coefficients = raw.coefficients,
                                      alphabest = alphabest,
                                      lambdabest = lambdabest)
  
  # Extracting elements of results
  objective <- results$objective
  xx <- results$xx
  raw.coefficients <- results$raw.coefficients
  coefficients <- results$coefficients
  a0 <- results$a0
  a00 <- results$a00
  raw.fitted.values <- results$raw.fitted.values
  raw.fitted.values.class <- results$raw.fitted.values.class
  fitted.values <- results$fitted.values
  fitted.values.class <- results$fitted.values.class
  
  ## OUTPUTS
  # Case: binomial
  if(family == "binomial"){
    output <- list(objective = objective, 
                   xx = xx, # TEMP REMOVE ME
                   best = sort(indexbest), 
                   raw.wt = raw.wt, 
                   wt = wgt, 
                   a00 = a00, 
                   raw.coefficients = raw.coefficients, 
                   a0 = a0, 
                   coefficients = coefficients, 
                   alpha = alphabest,
                   alphaw = alphaw_best,
                   lambda = lambdabest, 
                   lambdaw = lambdaw_best, 
                   num.nonzerocoef = num.nonzerocoef, 
                   h = h, 
                   raw.residuals = drop(raw.residuals), 
                   residuals = drop(reweighted.residuals), 
                   fitted.values = drop(fitted.values), 
                   fitted.values.class = drop(fitted.values.class), # NEW
                   raw.fitted.values = drop(raw.fitted.values), 
                   raw.fitted.values.class = drop(raw.fitted.values.class), # NEW
                   classnames = classnames, 
                   classsize = ntab, 
                   inputs = inputs, 
                   indexall = indexall, 
                   call = sys.call(),
                   alphas = alphas,  # NEW: Added the alphas that were used
                   lambdas = lambdas, # NEW: Added the lambdas that were used
                   ic_type = ic_type, # NEW: Added ic_type used
                   ic_type_reweighted = ic_type_reweighted) # NEW: Added ic_type_reweighted used
    
  } else if (family == "gaussian"){
    output <- list(objective = objective, 
                   best = sort(indexbest), 
                   raw.wt = raw.wt, 
                   wt = wgt, 
                   a00 = a00, 
                   raw.coefficients = raw.coefficients, 
                   a0 = a0, 
                   coefficients = coefficients, 
                   alpha = alphabest,
                   alphaw = alphaw_best,
                   lambda = lambdabest, 
                   lambdaw = lambdaw_best, 
                   num.nonzerocoef = num.nonzerocoef, 
                   h = h, 
                   raw.residuals = drop(raw.residuals), 
                   residuals = drop(reweighted.residuals), 
                   fitted.values = drop(fitted.values), 
                   raw.fitted.values = drop(raw.fitted.values), 
                   raw.rmse = raw.rmse, 
                   rmse = reweighted.rmse, 
                   inputs = inputs, 
                   indexall = indexall, 
                   call = sys.call(),
                   alphas = alphas,  # NEW: Added the alphas that were used
                   lambdas = lambdas)  # NEW: Added the lambdas that were used
  }
  class(output) <- "enetLTS"
  output$call <- matchedCall
  return(output)
}