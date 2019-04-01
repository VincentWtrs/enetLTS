enetLTS <- function(xx, yy, family = c("gaussian", "binomial"), alphas, 
                    lambdas, lambdaw, hsize = 0.75, intercept = TRUE, nsamp = 500, 
                    s1 = 10, nCsteps = 20, nfold = 5, seed = NULL, plot = TRUE, 
                    repl = 5, para = TRUE, ncores = 1, del = 0.0125, tol = -1e+06, 
                    scal = TRUE, type = c("response", "class"), ic_type = NULL, type_lambdaw = "min"){
  
  #### UPDATED VERSION ###
  # NEW: ic_type test test
  # NEW: type_lambdaw: choosing lambda.min or lambda.1se for the reweighting step.
  
  ## INPUT
  # xx: Predictor matrix (numeric, hence dummies need to be constructed beforehand)
  # yy: Outcome variable, for binomial factor CODED as 0/1!
  # family: The outcome variable distribution family (R Documentation states error distribution, this is wrong)
  # alphas: Sequence of alphas [0, 1] to fit the models, doesn't have to be in any order.
  # lambdas: Positive sequence of lambdas to fit the models, doesn't have to be in any order.
  # lamdbdaw: Positive sequence of lambdas for the reweighted fits
  # hsize: Fraction of data to be used for fitting. Default: 0.75, pick lower for more conservative fitting
  # intercept: Logical if intercept needs to be included. Best to use default (TRUE), and keep intercept out of xx!
  # scal: Logical if scaling of the h-subsamples needs to be performed, default: TRUE.
  # nsamp: The amount of elemental subsets (size = 3 or 4) to generate, default: 500.
  # s1: The amount of promising h (?) subsets to keep after running 2 C-steps on the nsamp elemental subsets.
  # nCsteps: A maximum amount of Csteps that the while-loops have to go through (to prevent stuck while-loops, default: 20).
  # nfold: Amount of folds to be used in the k-fold cross validation, default: 5.
  # seed: Initial seed for the random number generator for determining elemental subsets.
  # plot: Logical, if the alpha-lambda-Error tile plot needs to be made, default: TRUE.
  # repl: Amount of replications (repetitions) of the cross validation procedure, default: 5.
  # para: parallel fitting of the folds in the cross validation structure.
  # ncores: For parallel computation of certain parts (different than the CV structure!). NOT AVAILABLE ON WINDOWS.
  # del: 1-del: is the quantile to give 0-weight to an observation (hence this is flagged as outlier) in the reweighing step!
  # tol: Tolerance for stopping while loops in the C-step procedures, default: -1e+06
  # type: type of predictions required, default: c("response", "class")
  
  print("USING UPDATED VERSION OF enetLTS: WITH THE INTERCEPT FIX") # So I can see that the new function is effectively called
  ########################
  
  matchedCall <- match.call() # Use match.call() if there are a lot of optional arguments
  matchedCall[[1]] <- as.name("enetLTS")
  
  family <- match.arg(family)
  type <- match.arg(type)
  xx <- addColnames(as.matrix(xx)) # Adds column names
  nc = dim(yy) # Dimension of outcome matrix (n x 1)
  if (is.null(nc)) { # If yy is a vector, force it to matrix
    yy <- as.matrix(yy)
  }
  
  n <- nrow(xx) # Amount of observations
  p <- ncol(xx) # Amount of TRUE predictors (EXCLUDING INTERCEPT)
  h <- floor((n + 1) * hsize) # Size of outlier free set 
  
  ## Input parameters checks
  # CV repetitions (require positive)
  if (repl <= 0) {
    stop("repl has to be a positive number")
  }
  
  # Max C-step amount (require positive)
  if (nCsteps <= 0) {
    stop("nCsteps has to be a positive number")
  }
  
  # Output/Family consistency (class requires binomial)
  if (type == "class" & family == "gaussian") {
    stop("class type is not available for gaussian family")
  }
  
  # CPU cores (requires valid and <= than available cores)
  ncores <- rep(ncores, length.out = 1)
  if (is.na(ncores)) {
    ncores <- detectCores()
  }
  if (!is.numeric(ncores) || is.infinite(ncores) || ncores < 1) {
    ncores <- 1
    warning("invalid value of 'ncores'; using default value")
  } else {
    ncores <- as.integer(ncores)
  }
  
  # Information criterion (IC) (require nfold, repl to be 1)
  if(!is.null(ic_type)){
    print("Information criterion supplied: the cross-validation input parameters (nfold, repl) are ignored")
    nfold <- 1
    repl <- 1
  }
  
  ## Binomial response preparation
  if (family == "binomial") {
    rnames <- rownames(yy) # Saving rownames (TO DO: why?)
    rownames(yy) <- 1:nrow(yy) # Giving numbers as rownames
    y_actual = as.factor(yy) # Converting to factor
    ntab = table(y_actual) # Frequency table for y
    minclass = min(ntab) # Saving minority outcome class
    classnames = names(ntab) # Naming outcome class
  }
  
  ## alphas sequence
  # No given alphas: 41 evenly spaces values [0 ; 1]
  if (missing(alphas)) {
    alphas <- seq(from = 0, to = 1, length = 41)
  }
  
  # If given (user-supplied) alphas:
  alphas <- sort(alphas, decreasing = FALSE) # Increasing ordering
  
  wh <- (alphas < 0 | alphas > 1) # Checking if admissable
  if (sum(wh) > 0) { # Bit hacky but works
    stop("alphas can take the values only between 0 and 1")
  }
  
  ## lambdas sequence
  #################################
  # Missing and Gaussian: REMOVED #
  #################################
  
  
  if (missing(lambdas) & family == "binomial") { # Binomial case
    # Point-biserial method
    l00 <- lambda00(x = xx, 
                    y = yy, 
                    normalize = scal, 
                    intercept = intercept)
    lambdas <- seq(l00, 0, by = -0.025 * l00) # DECREASING SEQUENCE (HIGH REGULARIZATION FIRST) # TO DO: CHECK PROBABLY MORE LOGIC TO DO OTHER WAY AROUND TO GET OUTLIERS SINCE OTHERWISE THEY WILL BE HIDDEN IN THE SPAGHETTI (WARM START SEQUENCE) WHERE THE DIMENSIONALITY COLLAPSES AND THE OUTLIERS WILL NOT BE FOUND ANWAY...
    print("Sorting lambdas decreasingly: NEW FEATURE, seems more logical ")
    lambdas <- sort(lambdas, decreasing = FALSE) # NEW: Sorting INCREASING 
  }
  
  ## Data Preparation: robust centering and scaling
  sc <- prepara(x = xx, 
                y = yy, 
                family = family, # If "binomial": y not touched!
                robu = 1) # 
  # Extracting normalized (x,y)
  x <- sc$xnor # X is normalized (centering + scaling)
  y <- sc$ycen # In "binomial" case: no change!
  
  ## C-STEPS
  WarmCstepresults <- warmCsteps(x = x, 
                                 y = y, 
                                 h = h, 
                                 n = n, 
                                 p = p, 
                                 family = family, 
                                 alphas = alphas, # Calling with whole alphas sequence!
                                 lambdas = lambdas, # Calling with whole lambdas sequence!
                                 hsize = hsize, 
                                 nsamp = nsamp, 
                                 s1 = s1, 
                                 csteps = nCsteps, # Note differing argument name here (also note: small c in csteps and big one in nCsteps!)
                                 nfold = nfold, 
                                 para = para, 
                                 ncores = ncores, 
                                 tol = tol, 
                                 scal = scal, 
                                 seed = seed)
  indexall <- WarmCstepresults$indexall # Extracting indices of observations for all alpha-lambda combinations
  # NOTE: indexall is a 3-dim ARRAY (h x #lambdas x #alphas) (e.g. 75 x 41 x 20)
  # ... calling indexall[1, , ] will give the indices of the first observations for
  # ... all alpha-lambda values. Calling indexall[1, , 1] will give a vector of all
  # ... indices of the FIRST observations for all (e.g.) 41 lambda values (empty index)
  # ... for the first alpha value (third index in the indexall[] call)
  
  
  ## FINDING BEST ALPHA-LAMBDA VALUES
  # Case: 1x1 tuning grid (trivial)
  if ((length(alphas) == 1) & (length(lambdas) == 1)) {
    if(isTRUE(plot)){
      warning("There is no meaning to see plot for a single combination of lambda and alpha")
    }
    indexbest <- drop(indexall) # Dropping the matrix structure
    alphabest <- alphas # Obviously the only one is the best
    lambdabest <- lambdas # Obviously the only one is the best
  } 
  
  # Case: bigger tuning grid
  if ((length(alphas) > 1) | (length(lambdas) > 1)) {
    # NOTE: the results are called CV results can be from a IC search
    CVresults <- cv.enetLTS(index = indexall, 
                            xx = x, 
                            yy = y, 
                            family = family, 
                            h = h, 
                            alphas = alphas, 
                            lambdas = lambdas, 
                            nfold = nfold, 
                            repl = repl, 
                            ncores = ncores, 
                            plot = plot, 
                            ic_type = ic_type) # NEW: ic_type ## NEW(2): remove ic_type # NEW(3) Think we can remove this due to default NULL
    # Gathering results from search
    indexbest <- CVresults$indexbest
    alphabest <- CVresults$alphaopt
    lambdabest <- CVresults$lambdaopt
    evalCritCV <- CVresults$evalCrit
  }
  
  
  ### RESCALING STEP (Theoretically we should be in an outlier-free world now!)
  # Nonrobust scaling based on the outlier-free set
  if (isTRUE(scal)) { # TO DO: WHY IS THIS DONE?
    scl <- prepara(x = xx, 
                   y = yy, 
                   family = family, 
                   index = indexbest, # Indexes provided s.t. the normalization should happen on the outlier-free set
                   robu = 0) # Nonrobust because now in outlier-free world... (given by the indexbest index)
    xs <- scl$xnor # Normalized X
    ys <- scl$ycen # Centered y (not for binomial, confusing)
    
    # Fitting glmnet on best subset
    ## NEW: branch based on binomial
    if (family == "binomial") {
      fit <- glmnet(x = xs[indexbest, ],
                    y = ys[indexbest, ],
                    family = family,
                    alpha = alphabest,
                    lambda = lambdabest,
                    standardize = FALSE, # Because already done
                    intercept = TRUE) # Because in logit standardize and intercept don't interact in the same way 
    } else if (family == "gaussian") {
    fit <- glmnet(x = xs[indexbest, ], 
                  y = ys[indexbest, ], 
                  family = family, 
                  alpha = alphabest, 
                  lambda = lambdabest, 
                  standardize = FALSE, # Because already done
                  intercept = FALSE) 
  }
    
    ## Case: Binomial
    # Intercept handling
    if (family == "binomial") {
      if (isFALSE(intercept)) {
        a00 <- 0 # NEW: changed this from weird a00 <- if (intercept == FALSE) {0} style of notation
      } else {
        #a00 <- drop(fit$a0 - as.vector(as.matrix(fit$beta)) %*% (scl$mux / scl$sigx)) # NEW: same
        a00 <- fit$a0 # NEW: I think this is the way
      }
      
      # Extracting raw results
      raw.coefficients <- drop(as.matrix(fit$beta) / scl$sigx) # TO DO Does this actually hold for logit??
      beta_with_int <- cbind(a00, as.matrix(fit$beta) # NEW
      #raw.residuals <- -(ys * xs %*% as.matrix(fit$beta)) + log(1 + exp(xs %*% as.matrix(fit$beta))) # OLD
      raw.residuals <- -(ys * cbind(1, xs) %*% beta_with_in) + log(1 + exp(cbind(1, xs) %*% beta_with_int)) # NEW: cbind(1, xs) and beta_with_int
      raw.wt <- weight.binomial(x = xx, 
                                y = yy, 
                                beta = c(a00, raw.coefficients), 
                                intercept = intercept, 
                                del = del)
      
      # Normalizing using outlier-free data (raw.wt == 1) # TO DO: is raw.wt different from best index?
      sclw <- prepara(x = xx, 
                      y = yy, 
                      family = family, 
                      index = which(raw.wt == 1), # Only for those which are outlier-free
                      robu = 0)
      xss <- sclw$xnor 
      yss <- sclw$ycen
      
      # Tuning lambdaw (reweighted lambda)
      if (missing(lambdaw)) { # Case: no lambdaw given (DEFAULT)
        lambdaw_fit <- cv.glmnet(x = xss[which(raw.wt == 1), ], # NEW: changed name to lambdaw -> lambdaw_fit
                                 y = yss[which(raw.wt == 1)], 
                                 family = family, 
                                 nfolds = 5, 
                                 alpha = alphabest, 
                                 standardize = FALSE, 
                                 intercept = TRUE, # NEW: changed this to true
                                 type.measure = "deviance")
        # Note in case of no lambdaw given: it just uses the efficient algorithms!
        
        
        # Case: lambdaw given by user (unlikely)
      } else if (!missing(lambdaw) & length(lambdaw) == 1) { # Only single lambdaw given
        lambdaw <- lambdaw 
      } else if (!missing(lambdaw) & length(lambdaw) > 1) { # Multiple lambdaw given
        lambdaw_fit <- cv.glmnet(x = xss[which(raw.wt == 1), ], # NEW: changed name to lambdaw -> lambdaw_fit 
                                 y = yss[which(raw.wt == 1)], 
                                 family = family, 
                                 lambda = lambdaw, 
                                 nfolds = 5, 
                                 alpha = alphabest, 
                                 standardize = FALSE, 
                                 intercept = TRUE, # NEW: changed to this to true 
                                 type.measure = "deviance") # NEW: changed from "MSE" to "deviance" for the Gaussian case it is the same anyways
      }
      
      # Choosing lambda based on user-input (min vs. 1SE) # NEW
      if (type_lambdaw == "min") {
        lambdaw <- lambdaw_fit$lambda.min
      } else if (type_lambdaw == "1se") {
        lambdaw <- lambdaw_fit$lambda.1se
      }
      
      # Fitting using optimal lambdaw (reweighted!)
      fitw <- glmnet(x = xss[which(raw.wt == 1), ], 
                     y = yss[which(raw.wt == 1)], 
                     family = family, 
                     alpha = alphabest, 
                     lambda = lambdaw, 
                     standardize = FALSE, 
                     intercept = TRUE) # NEW/ changed this to true!
      
      # Intercept handling
      if (isFALSE(intercept)) {
        a0 <- 0 # NEW
      } else {
        #a0 <- drop(fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% (sclw$mux/sclw$sigx)) # NEW: REPLACED
        a0 <- fitw$a0
      }
      coefficients <- drop(as.matrix(fitw$beta)/sclw$sigx)
      wgt <- weight.binomial(x = xx, 
                             y = yy, 
                             beta = c(a0, coefficients), 
                             intercept = intercept, 
                             del = del)
      reweighted.residuals <- -(yy * cbind(1, xx) %*% c(a0, coefficients)) + log(1 + exp(cbind(1, xx) %*% c(a0, coefficients)))
    }
    
    ############################
    # Case: Gaussian: REmOVED #
    ############################
    
    
  } 
  
  ################################
  # Case: Scal == FALSE: REMOVED #
  ################################
  
  
  ## PREPARING OUTPUT
  # Counting number of nonzero coefficients
  num.nonzerocoef <- sum(coefficients != 0)
  
  # Adding intercept
  intercept <- isTRUE(intercept)
  if (intercept) { 
    xx <- addIntercept(x = xx)
  }
  if (intercept) {
    coefficients <- c(a0, coefficients)
    raw.coefficients <- c(a00, raw.coefficients)
  } else {
    coefficients <- coefficients
    raw.coefficients <- raw.coefficients
  }
  if (family == "binomial") {
    u <- xx %*% raw.coefficients # XBeta?
    raw.fitted.values <- if (type == "class") {
      ifelse(test = u <= 0.5, yes = 0, no = 1)
    } else if (type == "response"){
      1/(1 + exp(-u))
    }
    uu <- xx %*% coefficients
    fitted.values <- if (type == "class"){
      ifelse(test = uu <= 0.5, yes = 0, no = 1)
    } else if (type == "response"){
      1/(1 + exp(-uu))
    }
  }
  
  # TO DO: CHECK THIS
  if (family == "binomial"){
    objective <- h * (mean((-yy[indexbest] * (xx[indexbest, ] %*% coefficients)) + log(1 + exp(xx[indexbest, ] %*% coefficients))) + lambdabest * sum(1/2 * (1 - alphabest) * coefficients^2 + alphabest * abs(coefficients)))
  } else if (family == "gaussian"){
    objective <- h * ((1/2) * mean((yy[indexbest] - xx[indexbest, ] %*% coefficients)^2) + lambdabest * sum(1/2 * (1 - alphabest) * coefficients^2 + alphabest * abs(coefficients)))
  }
  if(intercept){
    coefficients <- coefficients[-1]
    raw.coefficients <- raw.coefficients[-1]
  } else {
    coefficients <- coefficients
    raw.coefficients <- raw.coefficients
  }
  inputs <- list(xx = xx, 
                 yy = yy, 
                 family = family, 
                 alphas = alphas, 
                 lambdas = lambdas, 
                 lambdaw = lambdaw, 
                 hsize = hsize, 
                 h = h, 
                 nsamp = nsamp, 
                 s1 = s1, 
                 nCsteps = nCsteps, 
                 nfold = nfold, 
                 intercept = intercept, 
                 repl = repl, 
                 para = para, 
                 ncores = ncores, 
                 del = del, 
                 scal = scal)
  
  ## OUTPUTS
  # Case: binomial
  if(family == "binomial"){
    output <- list(objective = objective, 
                   best = sort(indexbest), 
                   raw.wt = raw.wt, 
                   wt = wgt, 
                   a00 = a00, 
                   raw.coefficients = raw.coefficients, 
                   a0 = a0, 
                   coefficients = coefficients, 
                   alpha = alphabest, 
                   lambda = lambdabest, 
                   lambdaw = lambdaw, 
                   num.nonzerocoef = num.nonzerocoef, 
                   h = h, 
                   raw.residuals = drop(raw.residuals), 
                   residuals = drop(reweighted.residuals), 
                   fitted.values = drop(fitted.values), 
                   raw.fitted.values = drop(raw.fitted.values), 
                   classnames = classnames, 
                   classsize = ntab, 
                   inputs = inputs, 
                   indexall = indexall, 
                   call = sys.call())
    ############################  
    # Case: gaussian: REMOVED #
    ###########################
    
  }
  class(output) <- "enetLTS"
  output$call <- matchedCall
  output
}