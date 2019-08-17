enetLTS <- function(xx, yy, family = c("gaussian", "binomial"), alphas, 
                    lambdas, lambdaw, hsize = 0.75, intercept = TRUE, nsamp = 500, 
                    s1 = 10, nCsteps = 20, nfold = 5, seed = NULL, plot = TRUE, 
                    repl = 5, para = TRUE, ncores = 1, del = 0.0125, tol = -1e+06, 
                    scal = TRUE, type = c("response", "class"), ic_type = NULL, type_lambdaw = "min", 
                    ic_type_reweighted = NULL, simulation = FALSE) {
  
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
  
  # Output/Family consistency (class requires binomial) # TODO: TEMP: RE-ENABLE
  #if (type == "class" & family == "gaussian") {
  #  stop("class type is not available for gaussian family")
  #}
  
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
  
  # Check if in n < p problem such that at least some nonzero lambda is needed (only when no lambdas given)
  if(ncol(xx) > nrow(xx)) {
    min_lambda <- 1e-08
  } else {
    min_lambda <- 0 # 1e-08 #0 # TEMP NOT SETTING AT 0 BUT AT 1e-08
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
    lambdas <- seq(l00, min_lambda, by = -0.025 * l00) # DECREASING SEQUENCE (HIGH REGULARIZATION FIRST) # TO DO: CHECK PROBABLY MORE LOGIC TO DO OTHER WAY AROUND TO GET OUTLIERS SINCE OTHERWISE THEY WILL BE HIDDEN IN THE SPAGHETTI (WARM START SEQUENCE) WHERE THE DIMENSIONALITY COLLAPSES AND THE OUTLIERS WILL NOT BE FOUND ANWAY...
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
  if (is.null(ic_type)) {
    CVresults_list <- vector("list", length = 1) # TODO
  } else if (length(ic_type) >= 1){
    CVresults_list <- vector("list", length = length(ic_type)) # TODO
  }
  
  # Case: Tuning grid (more than single value of alpha and/or lambda)
  if ((length(alphas) > 1) | (length(lambdas) > 1)) {
    
    # Case: no ICs or single IC
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
      
      # Gathering results from the IC/CV search
      indexbest <- CVresults$indexbest
      alphabest <- CVresults$alphaopt
      lambdabest <- CVresults$lambdaopt
      evalCritCV <- CVresults$evalCrit
      
      # Case: multiple ICs provided
    } else if (length(ic_type) >= 1) {
      
      # Initialization of lists and indices
      #CVresults_list <- vector("list", length = length(ic_type)) # TODO CHECK WITH EARLIER INITIALIZATION ADDED
      names(CVresults_list) <- ic_type
      indexbest_list <- vector("list", length = length(ic_type))
      names(indexbest_list) <- ic_type
      print("structure of indexbest_list")
      print(str(indexbest_list))
      alphabest_list <- vector("list", length = length(ic_type))
      names(alphabest_list) <- ic_type
      lambdabest_list <- vector("list", length = length(ic_type))
      names(lambdabest_list) <- ic_type
      evalCritCV_list <- vector("list", length = length(ic_type))
      names(evalCritCV_list) <- ic_type
      
      # Looping over all required ICs
      for(i in 1:length(ic_type)) {
        ic_now <- ic_type[i] # Extract current IC
        print(paste0("Currently running the cv.enetLTS() function for: ", ic_type[i]))
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
        
        # Extracting results
        indexbest_list[[i]] <- CVresults_list[[i]]$indexbest
        alphabest_list[[i]] <- CVresults_list[[i]]$alphaopt
        lambdabest_list[[i]] <- CVresults_list[[i]]$lambdaopt
        evalCritCV_list[[i]] <- CVresults_list[[i]]$evalCrit
        print("I have finished a IC calculations and extracting pass") # TODO REMOVE
      }
    }
  }
  
  ### STEP: REWEIGHTING AND REFITTING
  ## Reweighting/Refitting with cross-validation
  if (is.null(ic_type)) {
    reweight_results <- enetLTS_reweight_results(xx = xx,
                                                 yy = yy,
                                                 family = family,
                                                 indexbest = indexbest,
                                                 alphabest = alphabest,
                                                 lambdabest = lambdabest,
                                                 h = h,
                                                 hsize = hsize,
                                                 nsamp = nsamp,
                                                 s1 = s1,
                                                 nCsteps = nCsteps,
                                                 nfold = nfold,
                                                 repl = repl,
                                                 para = para,
                                                 ncores = ncores,
                                                 ic_type = NULL, # TODO CHECK THIS IF WE CAN REMOVE THE THING
                                                 #lambdaw = lambdaw,
                                                 del = del,
                                                 intercept = intercept,
                                                 scal = scal,
                                                 type_lambdaw = type_lambdaw,
                                                 classnames = classnames,
                                                 ntab = ntab,
                                                 indexall = indexall,
                                                 alphas = alphas,
                                                 lambdas = lambdas)
  } else if (length(ic_type) >= 1) {
    
    # Initializing list to save results in
    reweight_results <- vector("list", length = length(ic_type)) # Initialization # TODO: MAKE MORE FLEXIBLE 
    
    # Looping over all ICs
    for (i in 1:length(ic_type)) {
      ic_type_now = ic_type[i] # Extracting current IC
      print("Entering the rweighting/results step for IC:")
      print(ic_type_now)
      
      reweight_results[[i]] <- enetLTS_reweight_results(xx = xx,
                                                        yy = yy,
                                                        family = family,
                                                        indexbest = indexbest_list[[i]],
                                                        alphabest = alphabest_list[[i]],
                                                        lambdabest = lambdabest_list[[i]],
                                                        h = h,
                                                        hsize = hsize,
                                                        nsamp = nsamp,
                                                        s1 = s1,
                                                        nCsteps = nCsteps,
                                                        nfold = nfold,
                                                        repl = repl,
                                                        para = para,
                                                        ncores = ncores,
                                                        ic_type = ic_type_now,
                                                        #lambdaw = lambdaw,
                                                        del = del,
                                                        intercept = intercept,
                                                        scal = scal,
                                                        type_lambdaw = type_lambdaw,
                                                        classnames = classnames,
                                                        ntab = ntab,
                                                        indexall = indexall,
                                                        alphas = alphas,
                                                        lambdas = lambdas)
    }
  }
  # OUTPUT
  return(reweight_results)
}