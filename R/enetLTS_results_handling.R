enetLTS_results_handling <- function(intercept, family, xx, a0, a00, coefficients, raw.coefficients, indexbest, alphabest, lambdabest) {
  ## Intercept handling
  intercept <- isTRUE(intercept)
  if (intercept) { 
    xx <- addIntercept(x = xx) # Adding a column of 1s
    coefficients <- c(a0, coefficients)
    raw.coefficients <- c(a00, raw.coefficients)
  } else if (!intercept) {
    coefficients <- coefficients
    raw.coefficients <- raw.coefficients
  }
  
  ### Fitted values BINOMIAL (NEW!)
  if (family == "binomial") {
    ## RAW FIT
    # Raw linear predictor (eta raw)
    u <- xx %*% raw.coefficients
    
    print(head(u, 100)) # TODO REMOVE
    
    # Raw fitted (predicted) probabilities
    raw.fitted.values <- 1/(1 + exp(-u))
    
    head(raw.fitted.values)
    
    # Raw fitted (predicted) classes based on 0.5 cutoff
    raw.fitted.values.class <- ifelse(test = raw.fitted.values > 0.5,
                                      yes = 1,
                                      no = 0)
    
    ## REWEIGHTED FIT
    # Reweighted linear predictor (eta reweighted)
    uu <- xx %*% coefficients
    
    # Reweighted fitted (predicted) probabilities
    fitted.values <- 1/(1 + exp(-uu))
    
    # Reweighted fitted (predicted) classes based on 0.5 cutoff
    fitted.values.class <- ifelse(test = fitted.values > 0.5,
                                  yes = 1,
                                  no = 0)
    ### Fitted values GAUSSIAN
  } else if (family == "gaussian") {
    raw.fitted.values <- xx %*% raw.coefficients
    fitted.values <- xx %*% coefficients
  }
  
  # TO DO: CHECK THIS
  if (family == "binomial") {
    objective <- h * (mean((-yy[indexbest] * (xx[indexbest, ] %*% coefficients)) + log(1 + exp(xx[indexbest, ] %*% coefficients))) + lambdabest * sum(1/2 * (1 - alphabest) * coefficients^2 + alphabest * abs(coefficients)))
  } else if (family == "gaussian") {
    objective <- h * ((1/2) * mean((yy[indexbest] - xx[indexbest, ] %*% coefficients)^2) + lambdabest * sum(1/2 * (1 - alphabest) * coefficients^2 + alphabest * abs(coefficients)))
  }
  if (intercept) {
    coefficients <- coefficients[-1] # Removing first element
    raw.coefficients <- raw.coefficients[-1]  # Removing first element
  } else {
    coefficients <- coefficients
    raw.coefficients <- raw.coefficients
  }
  if (family == "binomial") {
    output <- list(objective = objective,
                   xx = xx,
                   a00 = a00,
                   raw.coefficients = raw.coefficients,
                   a0 = a0,
                   coefficients = coefficients,
                   raw.fitted.values = raw.fitted.values,
                   raw.fitted.values.class = raw.fitted.values.class,
                   fitted.values = fitted.values,
                   fitted.values.class = fitted.values.class)
  } else if (family == "gaussian") {
    # OUTPUT TODO
    output <- list(objective = objective,
                   xx = xx,
                   a00 = a00,
                   raw.coefficients = raw.coefficients,
                   a0 = a0,
                   coefficients = coefficients,
                   raw.fitted.values = raw.fitted.values,
                   fitted.values = fitted.values)
  }
}