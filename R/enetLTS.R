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
    } 
    else if (length(ic_type) > 1) {
      print("I am in the multiple ic_type case! NOICE")
      # Inits + names
      CVresults_list <- vector("list", length = length(ic_type))
      names(CVresults_list) <- ic_type
      indexbest_list <- vector("list", length = length(ic_type))
      names(indexbest_list) <- ic_type
      alphabest_list <- vector("list", length = length(ic_type))
      names(alphabest_list) <- ic_type
      lambdabest_list <- vector("list", length = length(ic_type))
      names(lambdabest_list) <- ic_type
      evalCritCV <- vector("list", length = length(ic_type))
      names(evalCritCV) <- ic_type
      for(i in 1:length(ic_type)) {
        print(paste0("Currently running for: ", ic_type[i]))
        ic_now <- ic_type[i] # Extract IC
        CVresults_list[[i]] <- cv.enetLTS(index = indexall,
                                          xx = x,
                                          yy = y,
                                          family = family,
                                          alphas = alphas,
                                          lambdas = lambdas,
                                          nfold = nfold,
                                          repl = repl,
                                          ncores = ncores,
                                          plot = plot,
                                          ic_type = ic_now)
        #return(CVresults_list[[i]]) # TODO REMOVE ME
        print(CVresults_list[[i]]) # TEMP REMOVE ME TODO
        print('I am just before the indexbest extraction') # TODO REMOVE
        indexbest_list[[i]] <- CVresults_list[[i]]$indexbest
        print("I got past the indexbest extraction") # TODO REMOVE
        alphabest_list[[i]] <- CVresults_list[[i]]$alphaopt
        lambdabest_list[[i]] <- CVresults_list[[i]]$lambdaopt
        evalCritCV_list[[i]] <- CVresults_most[[i]]$evalCrit
        print("I have finished the first IC calculations")
      }
      
      if (simulation_run) {
        # TODO REMOVE
        print("I am just before the return statement now!")
        return(CVresults_list) # NEW TODO
      }
    }
  }
  
  #### FINAL FITTING OF THE MODEL USING THE OPTIMAL PARAMETERS (NON-REWEIGHTED FIT)
  ### Nonrobust scaling based on the outlier-free set
  ## Case: Scaling TRUE
  if (isTRUE(scal)) {
    
    # Scaling using the best index 'indexbest'
    scl <- prepara(x = xx, 
                   y = yy, 
                   family = family, 
                   index = indexbest, # Indexes provided s.t. the normalization should happen on the outlier-free set
                   robu = 0) # Nonrobust because now in outlier-free world... (given by the indexbest index)
    xs <- scl$xnor # Normalized X
    ys <- scl$ycen # Centered y (not for binomial, confusing)
    
    # Fitting: Case: Binomial (and scaling TRUE) 
    if (family == "binomial") {
      fit <- glmnet(x = xs[indexbest, ],
                    y = ys[indexbest, ],
                    family = "binomial",
                    alpha = alphabest,  # With the tuned alpha
                    lambda = lambdabest,  # ... and the tuned lambda
                    standardize = FALSE, # Because already done in the prepara case (we are providing y standardized and x standardized)
                    intercept = TRUE) # Because only for linear models the standardization makes the regression line go through the origin, not true for GLMs in general
      
      # Case: Gaussian (and scaling TRUE)
    } else if (family == "gaussian") {
      fit <- glmnet(x = xs[indexbest, ], 
                    y = ys[indexbest, ], 
                    family = "gaussian", 
                    alpha = alphabest, 
                    lambda = lambdabest, 
                    standardize = FALSE, # Because already done
                    intercept = FALSE)  # FALSE becasue for standardized data, a linear model goes through the origin
    }
    
    #### Results handling of final fit (non-reweighted)
    ### Intercept Handling
    ## Case: Binomial
    if (family == "binomial") {
      if (isFALSE(intercept)) {
        a00 <- 0 # NEW: changed this from weird a00 <- if (intercept == FALSE) {0} style of notation
      } else {
        #a00 <- drop(fit$a0 - as.vector(as.matrix(fit$beta)) %*% (scl$mux / scl$sigx)) # NEW: same # THIS WOULD MAKE SENSE IN LINEAR CASES
        a00 <- fit$a0 # NEW: I think this is the way (a0 is the way a glmnet object can have its intercept accessed)
        # NOTE: This only works if the model is fitted with a single lambda otherwise $a0 will give a vector of intercepts
      }
      
      # Extracting coefficients from final fit (non-reweighted)
      raw.coefficients <- drop(as.matrix(fit$beta) / scl$sigx) # This holds for all GLMs, invariance principle
      beta_with_int <- drop(as.matrix((coef(fit, s = lambdabest)))) # Adding s = lambdabest as redundancy, if no lambdabest would be supplied to glmnet, it would be able to display it for all lambdas
      # Note calling coef on a glmnet object will give those dcg sparse matrices. By forcing to matrix and dropping unnecessary dimensions we get a named vector
      #raw.residuals <- -(ys * xs %*% as.matrix(fit$beta)) + log(1 + exp(xs %*% as.matrix(fit$beta))) # OLD
      
      raw.residuals <- -(ys * cbind(1, xs) %*% beta_with_int) + log(1 + exp(cbind(1, xs) %*% beta_with_int)) # NEW: cbind(1, xs) and beta_with_int
      # NOTE: This should be the formula of the residuals
      raw.wt <- weight.binomial(x = xx, 
                                y = yy, 
                                beta = c(a00, raw.coefficients),  # Or beta_with_int
                                intercept = intercept, 
                                del = del)
      
      ### REWEIGHTING STEP
      # Again: Normalizing using outlier-free data (raw.wt == 1) # TODO: is raw.wt different from best index?
      sclw <- prepara(x = xx, 
                      y = yy, 
                      family = family, 
                      index = which(raw.wt == 1), # Only for those which are outlier-free
                      robu = 0)
      xss <- sclw$xnor 
      yss <- sclw$ycen
      
      ### REWEIGHTING STEP
      ## Tuning lambdaw
      # Case: no given lambdaw (Default)
      if (missing(lambdaw)) {
        print("We are here, where I will do the cva.glmnet")
        #reweighted_cv <- cv.glmnet(x = xss[which(raw.wt == 1), ], # NEW: changed name to lambdaw -> lambdaw_fit
        #                         y = yss[which(raw.wt == 1)], 
        #                         family = family, 
        #                         nfolds = 5, 
        #                         alpha = alphabest, 
        #                         standardize = FALSE, 
        #                         intercept = TRUE, # NEW: changed this to true
        #                         type.measure = "deviance")
        reweighted_cv <- cva.glmnet(x = xss[which(raw.wt == 1), ],
                                    y = yss[which(raw.wt == 1)],
                                    family = family,
                                    nfolds = 5,
                                    standardize = FALSE,
                                    intercept = TRUE)
        
        # Extracting the models best alpha index
        alphaw_index <- which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm)))
        
        # Assigning as the new alphabest
        alphabest <- reweighted_cv$alpha[alphaw_index]
        
        # Note in case of no lambdaw given: it just uses the efficient algorithms!
        
        # Maybe extract the used lambdaw also
        
        # NEW: ADDING IC-BASED REWEIGHTED TUNING # TO DO: fix this, we need to do the IC tuning at the end, get proper loss, etc (which is the loss now to take into account...)
        #if (!is.null(ic_type_reweighted)) {
        #  for(1:length(lambdas)){
        #    lambda
        #  }
        #  lambdaw_fit <- glmnet(x = xss[which(raw.wt == 1), ],
        #                        y = yss[which(raw.wt == 1)],
        #                        family = family,
        #                        alpha = alphabest,
        #                        )
        #}
        
        
        # Case: SINGLE lambdaw given by user (unlikely)
      } else if (!missing(lambdaw) & length(lambdaw) == 1) { # Only single lambdaw given
        lambdaw <- lambdaw
        
        # Case: Multiple lambdaws given by user (less unlikely)
      } else if (!missing(lambdaw) & length(lambdaw) > 1) { # Multiple lambdaw given
        reweighted_cv <- cva.glmnet(x = xss[which(raw.wt == 1), ], # NEW: changed name to lambdaw -> lambdaw_fit 
                                    y = yss[which(raw.wt == 1)], 
                                    family = family, 
                                    lambda = lambdaw, 
                                    nfolds = 5, 
                                    #alpha = alphabest, # TODO: Check this
                                    standardize = FALSE, 
                                    intercept = TRUE, # NEW: changed to this to true 
                                    type.measure = "deviance") # NEW: changed from "MSE" to "deviance" for the Gaussian case it is the same anyways
      }
      
      # Choosing lambda based on user-input (min vs. 1SE) # NEW
      if (type_lambdaw == "min") {
        print("Extracting the optimal lambdaw from cva.glmnet") # TODO: REMOVE
        #lambdaw <- lambdaw_fit$lambda.min
        #lambdaw <- reweighted_cv$lambda[which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm)))]  # TODO (TEMP) Correct, this is to see if it works 
        lambdaw <- reweighted_cv$modlist[[alphaw_index]]$lambda.min # NEW!
      } else if (type_lambdaw == "1se") {
        #lambdaw <- lambdaw_fit$lambda.1se
        #lambdaw <- reweighted_cv$lambda[which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm)))]  # TODO (TEMP) Correct, this is to see if it works
        lambdaw <- reweighted_cv$modlist[[alphaw_index]]$lambda.1se
      }
      
      # Fitting using optimal lambdaw (reweighted!)
      print("alphabest:") # TODO REMOVE
      print(alphabest) # TODO REMVOE
      print("lambdaw:") # TODO REMOVE
      print(lambdaw) # TODO REMOVE
      fitw <- glmnet(x = xss[which(raw.wt == 1), ], 
                     y = yss[which(raw.wt == 1)], 
                     family = family, 
                     alpha = alphabest, 
                     lambda = lambdaw, 
                     standardize = FALSE, 
                     intercept = TRUE) # NEW: changed this to true!
      # TODO REMOVE
      print("I am past fitting the fitw object")
      
      # Intercept handling
      if (isFALSE(intercept)) {
        a0 <- 0 # NEW
      } else {
        #a0 <- drop(fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% (sclw$mux/sclw$sigx)) # NEW: REPLACED
        a0 <- fitw$a0
      }
      
      coefficients <- drop(as.matrix(fitw$beta)/sclw$sigx)
      print("Structure of coefficients:") # TODO REMOVE
      print(str(coefficients)) # # TODO REMOVE
      print("Structure of a0:") # # TODO REMOVE
      print(str(a0)) # TODO REMOVE
      
      
      print("Calculated coefs")  # TODO REMOVE
      wgt <- weight.binomial(x = xx, 
                             y = yy, 
                             beta = c(a0, coefficients), 
                             intercept = intercept, 
                             del = del)
      print("ran weight binomial here")  # TODO REMOVE
      
      reweighted.residuals <- -(yy * cbind(1, xx) %*% c(a0, coefficients)) + log(1 + exp(cbind(1, xx) %*% c(a0, coefficients)))
      print("Calculated reweighed resids")  # TODO REMoVE
      
      # Case: Gaussian
    } else if (family == "gaussian") {
      a00 <- if (intercept == FALSE) 
        0
      else drop(scl$muy + fit$a0 - as.vector(as.matrix(fit$beta)) %*% (scl$mux/scl$sigx))
      
      raw.coefficients <- drop(as.matrix(fit$beta)/scl$sigx)
      raw.residuals <- yy - cbind(1, xx) %*% c(a00, raw.coefficients)
      raw.rmse <- sqrt(mean(raw.residuals^2))
      raw.wt <- weight.gaussian(resi = raw.residuals, 
                                ind = indexbest, 
                                del = del)$we
      sclw <- prepara(x = xx, 
                      y = yy, 
                      family = family, 
                      index = which(raw.wt == 1), 
                      robu = 0)
      xss <- sclw$xnor
      yss <- sclw$ycen
      
      if ((missing(lambdaw))) {
        reweighted_cv <- cv.glmnet(x = xss[which(raw.wt == 1), ],
                                   y = yss[which(raw.wt == 1)], 
                                   family = family, 
                                   nfolds = 5, 
                                   alpha = alphabest, 
                                   standardize = FALSE, 
                                   intercept = FALSE, 
                                   type.measure = "mse") # NEW: REMOVED $lambda.min here because we extract it later anyways
      }
      else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw <- lambdaw
      }
      else if (!missing(lambdaw) & length(lambdaw) > 1) {
        reweighted_cv <- cv.glmnet(x = xss[which(raw.wt == 1), ],
                                   y = yss[which(raw.wt == 1)], 
                                   family = family, 
                                   lambda = lambdaw, 
                                   nfolds = 5, 
                                   alpha = alphabest, 
                                   standardize = FALSE, 
                                   intercept = FALSE, 
                                   type.measure = "mse")
      }
      
      # NEW: Choosing lambda
      if (type_lambdaw == "min") {
        lambdaw <- lambdaw_fit$lambda.min 
      } else if (type_lambdaw == "1se") {
        lambdaw <- lambdaw_fit$lambda.1se
      }
      
      # Fitting using optimal lambdaw
      fitw <- glmnet(x = xss[which(raw.wt == 1), ], 
                     y = yss[which(raw.wt == 1)], 
                     family = "gaussian", 
                     alpha = alphabest, 
                     lambda = lambdaw, 
                     standardize = FALSE, 
                     intercept = FALSE)
      a0 <- if (intercept == FALSE) 
        0
      else drop(sclw$muy + fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% (sclw$mux/sclw$sigx))
      
      coefficients <- drop(as.matrix(fitw$beta)/sclw$sigx)
      reweighted.residuals <- yy - cbind(1, xx) %*% c(a0, coefficients)
      reweighted.rmse <- sqrt(mean(reweighted.residuals^2))
      wgt <- weight.gaussian(resi = reweighted.residuals, 
                             ind = raw.wt == 1, 
                             del = del)$we
    } # End of Case: Gaussian
    # End of isTRUE(scal)
  } else { # Think (!) this amounts to having isFALSE(scal):
    # Case: No scaling (unlikely)
    fit <- glmnet(x = x[indexbest, ], 
                  y = y[indexbest, ], 
                  family = family, 
                  alpha = alphabest, 
                  lambda = lambdabest, 
                  standardize = FALSE, 
                  intercept = FALSE)
    if (family == "binomial") {
      a00 <- if (intercept == FALSE) 
        0
      else drop(fit$a0 - as.vector(as.matrix(fit$beta)) %*% (sc$mux/sc$sigx))
      raw.coefficients <- drop(as.matrix(fit$beta)/sc$sigx)
      raw.residuals <- -(y * x %*% as.matrix(fit$beta)) + log(1 + exp(x %*% as.matrix(fit$beta)))
      raw.wt <- weight.binomial(x = xx, 
                                y = yy, 
                                beta = c(a00, raw.coefficients), 
                                intercept = intercept, 
                                del = del)
      # This is probably the default:
      if (missing(lambdaw)) {
        #reweighted_cv <- cv.glmnet(x = x[which(raw.wt == 1), ], 
        #                          y = y[which(raw.wt == 1)], 
        #                          family = family, 
        #                          nfolds = 10, 
        #                          alpha = alphabest, 
        #                          standardize = FALSE, 
        #                          intercept = FALSE, 
        #                          type.measure = "deviance") #$lambda.min # removed check check NEW
        reweighted_cv <- cva.glmnet(x = x[which(raw.wt == 1), ], 
                                    y = y[which(raw.wt == 1)], 
                                    family = family, 
                                    nfolds = 10, 
                                    #alpha = alphabest, # TODO (Check)
                                    standardize = FALSE, 
                                    intercept = FALSE)
        alphaw <- reweighted_cv$alpha[which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm)))]
        print(paste0("The reweighted alpha is: ", alphaw))
        
        
        
        
      }
      else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw <- lambdaw
      }
      # In case of multiple user supplied lambdaws
      else if (!missing(lambdaw) & length(lambdaw) > 1) {
        # NEW: saved model first, then extract the lambda we want
        reweighted_cv <- cv.glmnet(x = x[which(raw.wt == 1), ], 
                                   y = y[which(raw.wt == 1)], 
                                   family = family,
                                   lambda = lambdaw, 
                                   nfolds = 10, 
                                   alpha = alphabest, 
                                   standardize = FALSE, 
                                   intercept = FALSE, 
                                   type.measure = "deviance") # NEW: changed type.measure = "mse" to type.meaure = "deviance" 
        # NEW: Plotting the cross-validation curves
        
        
        # NEW: removed cv.glmnet(...)$lambda.min and given choice!
        if(type_lambdaw == "min") {
          lambdaw <- reweighted_fit$lambda.min
        } else if(type_lambdaw == "1se") {
          lambdaw <- reweighted_fit$lambda.1se
        }
      }
      
      # Fitting using optimal lambdaw
      fitw <- glmnet(x = x[which(raw.wt == 1), ], 
                     y = y[which(raw.wt == 1)], 
                     family = family, 
                     alpha = alphabest, 
                     lambda = lambdaw, 
                     standardize = FALSE, 
                     intercept = FALSE)
      a0 <- if(intercept == FALSE) 
        0
      else drop(fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% (sc$mux/sc$sigx))
      coefficients <- drop(as.matrix(fitw$beta)/sc$sigx)
      wgt <- weight.binomial(x = xx, 
                             y = yy, 
                             beta = c(a0, coefficients), 
                             intercept = intercept, 
                             del = del)
      reweighted.residuals <- -(yy * cbind(1, xx) %*% c(a0, coefficients)) + log(1 + exp(cbind(1, xx) %*% c(a0, coefficients)))
    } else if (family == "gaussian") { # Case: Gaussian (and no scaling)
      a00 <- if (intercept == FALSE) 
        0
      else drop(sc$muy + fit$a0 - as.vector(as.matrix(fit$beta)) %*% (sc$mux/sc$sigx))
      raw.coefficients <- drop(as.matrix(fit$beta)/sc$sigx)
      raw.residuals <- yy - cbind(1, xx) %*% c(a00, raw.coefficients)
      raw.rmse <- sqrt(mean(raw.residuals^2))
      raw.wt <- weight.gaussian(resi = raw.residuals, 
                                ind = indexbest, 
                                del = del)$we
      
      if (missing(lambdaw)) {
        reweighted_cv <- cv.glmnet(x = x[which(raw.wt == 1), ], 
                                   y = y[which(raw.wt == 1)], 
                                   family = family, 
                                   nfolds = 5, 
                                   alpha = alphabest, 
                                   standardize = FALSE, 
                                   intercept = FALSE, 
                                   type.measure = "mse") # $lambda.min # NEW REMOVED
        
      }
      else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw <- lambdaw
      }
      else if (!missing(lambdaw) & length(lambdaw) > 1) {
        reweighted_cv <- cv.glmnet(x = x[which(raw.wt == 1), ], 
                                   y = y[which(raw.wt == 1)], 
                                   family = family, 
                                   lambda = lambdaw, 
                                   nfolds = 5, 
                                   alpha = alphabest, 
                                   standardize = FALSE, 
                                   intercept = FALSE, 
                                   type.measure = "mse")
      }
      
      # NEW: lambda choosing part
      if (type_lambdaw == "min") {
        #lambdaw <- reweighted_cv$lambda.min
        lambdaw <- reweighted_cv$lambda[which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm)))]  # TODO (TEMP) Correct, this is to see if it works
        
        
      } else if (type_lambdaw == "1se") {
        #lambdaw <- reweighted_cv$lambda.1se
        lambdaw <- reweighted_cv$lambda[which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm)))]  # TODO (TEMP) correct this, this is to see if it works
        
      }
      
      # Fitting using optimal lambda
      fitw <- glmnet(x = x[which(raw.wt == 1), ], 
                     y = y[which(raw.wt == 1)], 
                     family = family, 
                     alpha = alphabest, 
                     lambda = lambdaw, 
                     standardize = FALSE, 
                     intercept = FALSE)
      a0 <- if(intercept == FALSE) 
        0
      else drop(sc$muy + fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% (sc$mux/sc$sigx))
      coefficients <- drop(as.matrix(fitw$beta)/sc$sigx)
      reweighted.residuals <- yy - cbind(1, xx) %*% c(a0, coefficients)
      reweighted.rmse <- sqrt(mean(reweighted.residuals^2))
      wgt <- weight.gaussian(resi = reweighted.residuals, 
                             ind = raw.wt == 1, 
                             del = del)$we
    }
  } # End scal == FALSE
  
  # Plotting CV Plot 
  plot(reweighted_cv)
  print(paste0("The optimal reweighted lambda is: ", lambdaw)) # Printing chosen lambdaw
  print(paste0("The optimal alpha used in the reweighting step is: ", alphabest)) # Printing chosen alpha
  # Note: I put this here at the end because otherwise I would need to repeat the statement multiple times
  
  ## PREPARING OUTPUT
  # Counting number of nonzero coefficients
  num.nonzerocoef <- sum(coefficients != 0)
  
  # Adding intercept
  intercept <- isTRUE(intercept)
  if (intercept) { 
    xx <- addIntercept(x = xx)
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
    fitted.values <- if (type == "class") {
      ifelse(test = uu <= 0.5, yes = 0, no = 1)
    } else if (type == "response") {
      1/(1 + exp(-uu))
    }
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
                   call = sys.call(),
                   alphas = alphas,  # NEW: Added the alphas that were used
                   lambdas = lambdas)  # NEW: Added the lambdas that were used
    
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
                   lambda = lambdabest, 
                   lambdaw = lambdaw, 
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
  output
}