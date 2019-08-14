enetLTS <- function(xx, yy, family = c("gaussian", "binomial"), alphas, 
                    lambdas, lambdaw, hsize = 0.75, intercept = TRUE, nsamp = 500, 
                    s1 = 10, nCsteps = 20, nfold = 5, seed = NULL, plot = TRUE, 
                    repl = 5, para = TRUE, ncores = 1, del = 0.0125, tol = -1e+06, 
                    scal = TRUE, type = c("response", "class"), ic_type = NULL, type_lambdaw = "min", ic_type_reweighted = NULL, simulation = FALSE){
  
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
  
  # Loading libraries
  library(glmnetUtils) # For cva.glmnet(): Tuning of alpha and lambda
  
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
  # Case: No user given lambdas and gaussian
  if (missing(lambdas) & family == "gaussian") {
    # Bivariate winsorization method
    l0 <- robustHD::lambda0(xx, yy, normalize = scal, intercept = intercept)
    lambdas <- seq(l0, 0, by = -0.025 * l0) # DECREASING SEQUENCE (HIGHEST REGULARIZATION FIRST) # TO DO: CHECK PROBABLY MORE LOGIC TO DO OTHER WAY AROUND TO GET OUTLIERS SINCE OTHERWISE THEY WILL BE HIDDEN IN THE SPAGHETTI (WARM START SEQUENCE) WHERE THE DIMENSIONALITY COLLAPSES AND THE OUTLIERS WILL NOT BE FOUND ANWAY...
    
    # Case: No user given lambdas and binomial
  } else if (missing(lambdas) & family == "binomial") { # Binomial case
    # Point-biserial method
    l00 <- lambda00(x = xx, 
                    y = yy, 
                    normalize = scal, 
                    intercept = intercept)
    lambdas <- seq(l00, 0, by = -0.025 * l00) # DECREASING SEQUENCE (HIGH REGULARIZATION FIRST) # TO DO: CHECK PROBABLY MORE LOGIC TO DO OTHER WAY AROUND TO GET OUTLIERS SINCE OTHERWISE THEY WILL BE HIDDEN IN THE SPAGHETTI (WARM START SEQUENCE) WHERE THE DIMENSIONALITY COLLAPSES AND THE OUTLIERS WILL NOT BE FOUND ANWAY...
    #lambdas <- sort(lambdas, decreasing = FALSE) # NEW: Sorting INCREASING 
    # TODO (Seeing if keeping in original order fixes the problems)
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
  
  # INITIALIZING # TODO TEMP REMOVE ME / IMPROVE ME
  if(is.null(ic_type)){
    CVresults_list <- vector("list", length = 1) # TODO
  } else if (length(ic_type) >= 1){
    CVresults_list <- vector("list", length = length(ic_type)) # TODO
  }
  
  # Case: bigger tuning grid
  if ((length(alphas) > 1) | (length(lambdas) > 1)) {
    if ((is.null(ic_type)) | (length(ic_type) == 1)) {
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
    }  else if (length(ic_type) > 1) {
      print("I am in the multiple ic_type case! NOICE")
      
      # Initialization of lists and indices
      #CVresults_list <- vector("list", length = length(ic_type)) # TODO CHECK WITH EARLIER INITIALIZATION ADDED
      names(CVresults_list) <- ic_type
      indexbest_list <- vector("list", length = length(ic_type))
      names(indexbest_list) <- ic_type
      alphabest_list <- vector("list", length = length(ic_type))
      names(alphabest_list) <- ic_type
      lambdabest_list <- vector("list", length = length(ic_type))
      names(lambdabest_list) <- ic_type
      evalCritCV_list <- vector("list", length = length(ic_type))
      names(evalCritCV_list) <- ic_type
      
      # Looping over all requires ICs
      for(i in 1:length(ic_type)) {
        ic_now <- ic_type[i] # Extract current IC
        print(paste0("Currently running for: ", ic_type[i]))
        CVresults_list[[i]] <- cv.enetLTS(index = indexall,
                                          xx = x, # x Is the normalized data, xx is raw
                                          yy = y, # y Is the normalized data, yy is raw
                                          family = family,
                                          alphas = alphas,
                                          lambdas = lambdas,
                                          nfold = nfold,
                                          repl = repl,
                                          ncores = ncores,
                                          plot = plot,
                                          ic_type = ic_now)
        #return(CVresults_list[[i]]) # TODO REMOVE ME
        #print(CVresults_list[[i]]) # TEMP REMOVE ME TODO
        #print('I am just before the indexbest extraction') # TODO REMOVE
        
        # Extracting results
        indexbest_list[[i]] <- CVresults_list[[i]]$indexbest
        #print("I got past the indexbest extraction") # TODO REMOVE
        alphabest_list[[i]] <- CVresults_list[[i]]$alphaopt
        print(paste0("I am printing the extracted best alpha: ", alphabest_list[[i]])) # TODO REMOVE
        lambdabest_list[[i]] <- CVresults_list[[i]]$lambdaopt
        evalCritCV_list[[i]] <- CVresults_list[[i]]$evalCrit
        
        print("I have finished a IC calculations and extracting pass")
      }
      
      if (simulation) {
        # TODO REMOVE
        print("I am just before the return statement now!")
        #return(CVresults_list) # CHECK WHAT'S IN THERE
      }
    }
  }
  
  ## SIMULATION-RELATED CODE
  if(simulation){
    output_list <- vector("list", length = length(ic_type)) # Used to loop the results later
  } else {
    output_list <- vector("list", length = 1) # A single length list that will later be extracted so the original output is as before again
  }
  
  
  for (i in 1:length(ic_type)){
    if (simulation) { # NEW # TODO CHECK / SIMPLIFY THIS
      # Results extraction from list
      print("I am extracting the results from the CVresults_list")
      indexbest = indexbest_list[[i]]
      alphabest <- alphabest_list[[i]]
      lambdabest <- lambdabest_list[[i]]
      evalCritCV <- evalCritCV_list[[i]]
    }
    
    ### STEP: REWEIGHTING AND REFITTING
    ## Running function
    # Case for missing lambdaw (call function without)
    if(missing(lambdaw)){
      print('Printing lambdabest before calling enetLTS_reweighting_refitting')
      print(lambdabest)
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
    
    print("I am past running the enetLTS_reweighting_refitting function") # TODO REMOVE 
    
    ## Extracting results
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
    
    print("printing raw coefs just after extraction")
    print(raw.coefficients)
    
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
    
    ## Intercept handling
    intercept <- isTRUE(intercept)
    if (intercept) { 
      #xx <- addIntercept(x = xx) # Adding a column of 1s # TEMPO REMOVE THIS BECAUSE I WANT TO TEST SOMETHING
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
      
      print("structure of xx")
      print(str(xx))
      
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
    
    inputs <- list(xx = xx, 
                   yy = yy, 
                   family = family, 
                   alphas = alphas, 
                   lambdas = lambdas, 
                   #lambdaw = lambdaw, # TODO REMOVED BECAUSE CAUSING SOME ISSUES IN NEW FUNCTION
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
    output[[i]] <- output
    print("I have done everything for the IC: ")
    print(ic_type[i])
  }
  
  if(!simulation) {
    output <- output[[1]] # In this case: just single run, hence just give it the single element
  }
}