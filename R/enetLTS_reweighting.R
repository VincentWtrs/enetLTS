enetLTS_reweighting_refitting <- function (xx, yy, family, indexbest, alphabest, lambdabest, lambdaw, intercept) {
  
  # In this part we will get a final fit. Common practice in machine learning, after the optimal hyperparameters have been found
  # ... to take this hyperparameter and fit it to all data. Here, this case, it is slightly more complex: there is (i) the tuned model
  # ... that is "non-reweighted", for which a final fit can be obtained. Secondly (ii) there is the reweighted tuned model for which
  # ... a final fit can be obtained as well. The issue is that not all data is used but only the 'supposed' outlier-free set, of course.
  # ... The reweighting step also complicates things even further since the amount of outliers will be reduced.
  
  ## INPUTS
  # xx: data
  # yy: data
  # family: "binomial", "gaussian" is WiP
  # indexbest: the found outlier-free indices
  # alphabest: Best raw alpha
  # lambdabest: Best raw lambda
  # intercept: If intercept to be used or not (FALSE is WiP)
  
  #### FINAL FITTING OF THE MODEL USING THE OPTIMAL PARAMETERS (NON-REWEIGHTED FIT)
  
  ### STEPS FOR SCAL: TRUE
  ## STEP 1. Scaling (if scal TRUE) using the best indices found (related to best hyperparameter), hence nonrobust (because outlier-free world)
  ## STEP 2. Fitting of ordinary elastic net (glmnet) on all non-outlying data with best hyperparams and scaled outlier-free dataset
  ## STEP 3. Handling of the results from the final (non-reweighted) fit.
  # STEP 3a. Intercept handling of (non-reweighted) fit
  # STEP 3b. Coefficients handling of (non-reweighred) fit
  # STEP 3c. Residuals handling of (non-reweighted) fit
  #---------------- REWEIGHTING ----------------------#
  ## STEP 4. Calculating raw weights (raw.wt)
  ## STEP 5. Scaling (if scal TRUE) using the best REWEIGHTED indices found found (related to best hyperparamters), hence nonrobust but potentially less conservative
  ## STEP 6. HYPERPARAMETER TUNING
  # STEP 6a. Fitting all hyperparameters their models (i.e. calling cva.glmnet())
  # STEP 6b. CHOOSING OPTIMAL ALPHA (JUST MIN)
  # STEP 6c. CHOOSING OPTIMAL LAMBDA
  ## STEP 7
  ## STEP 8b. INTERCEPT HANDLING OF REWEIGHTED FIT
  ## STEP 8b. COEFFICIENT HANDLING OF REWEIGHTED FIT
  
  ### STEP 0. Printing some warnings
  if (family == "gaussian") {
    warning("Currently the functionality is not properly tested for the Gaussian family, this might throw errors or not work as desired")
  }
  if (isFALSE(intercept)) {
    warning("Currently the functionality is not properly tested for the no-intercept situation (FALSE), this might throw errors or not work as expected")
  }
  
  ### STEP 1. Scaling of data
  if (isTRUE(scal)) {
    scl <- prepara(x = xx, 
                   y = yy, 
                   family = family, 
                   index = indexbest, # Indices provided s.t. the normalization should happen on the outlier-free set (=='outlier-free world')
                   robu = 0) # Nonrobust because now in outlier-free world, as provided by the indices above
    xs <- scl$xnor # Normalized X
    ys <- scl$ycen # Centered y (normally, but since family is passed on, no centering is done for Binomial!)
    
    ### STEP 2: Final fitting elastic net with optimal hyperparameters and outlier-free set (non-reweighted)
    ## STEP 2: Final fitting: Case Binomial
    if (family == "binomial") {
      fit <- glmnet(x = xs[indexbest, ],
                    y = ys[indexbest, ],
                    family = "binomial",
                    alpha = alphabest,  # With the tuned alpha
                    lambda = lambdabest,  # ... and the tuned lambda
                    standardize = FALSE, # Because already done in the prepara case (we are providing y standardized and x standardized)
                    intercept = TRUE) # Because only for linear models the standardization makes the regression line go through the origin, not true for GLMs in general
      
      ## STEP 2: Final fitting: Case Gaussian
    } else if (family == "gaussian") {
      fit <- glmnet(x = xs[indexbest, ], 
                    y = ys[indexbest, ], 
                    family = "gaussian", 
                    alpha = alphabest, 
                    lambda = lambdabest, 
                    standardize = FALSE, # Because already done
                    intercept = FALSE)  # FALSE becasue for standardized data, a linear model goes through the origin!
    }
    
    ### STEP 3: HANDLING RESULTS OF THE FINAL (NON-REWEIGHTED) FIT
    ## STEP 3a. Intercept Handling
    # STEP 3a. Intercept handling: Case Binomial
    if (family == "binomial") { # THIS GOES ON FOR A LONG WHILE
      if (isFALSE(intercept)) {
        a00 <- 0 # NEW: Changed this from weird a00 <- if (intercept == FALSE) {0} style of notation, this is clearer, does the same
      } else if (isTRUE(intercept)) {
        a00 <- fit$a0 # NEW: I think this is the way (a0 is the way a glmnet object can have its intercept accessed)
        #a00 <- drop(fit$a0 - as.vector(as.matrix(fit$beta)) %*% (scl$mux / scl$sigx)) # NEW: REMOVED SINCE THIS WOULD ONLY! MAKE SENSE IN LINEAR CASES
        # NOTE: This only works if the model is fitted with a single lambda otherwise $a0 will give a vector of intercepts
      }
      
      # STEP 3b. Coefficients handling of (non-reweighted) fit
      raw.coefficients <- drop(as.matrix(fit$beta) / scl$sigx) # Getting original coefs (w.r.t. unstandardized data)
      # NOTE: Formula above holds for all GLMs due to invariance principle. Note that this must be done since the data provided was standardized before the function! I think here we could just run glmnet as-is as well...
      beta_with_int <- drop(as.matrix((coef(fit, s = lambdabest)))) # Adding s = lambdabest as redundancy, if no lambdabest would be supplied to glmnet, it would be able to display it for all lambdas
      # Note calling coef on a glmnet object will give those dcg sparse matrices. By forcing to matrix and dropping unnecessary dimensions we get a named vector
      
      # STEP 3c. Residuals handling of (non-reweighted) fit
      raw.residuals <- -(ys * cbind(1, xs) %*% beta_with_int) + log(1 + exp(cbind(1, xs) %*% beta_with_int)) # NEW: cbind(1, xs) AND beta_with_int TO TAKE INTERCEPT INTO ACCOUNT
      #raw.residuals <- -(ys * xs %*% as.matrix(fit$beta)) + log(1 + exp(xs %*% as.matrix(fit$beta))) # OLD FORMULA, NOT CORRECT FOR BINOMIAL SINCE IT DOESN'T TAKE INTERCEPT INTO ACCOUNT
      
      ##################### REWEIGHTING STEP #####################
      
      ## STEP 4: CALCULATING RAW WEIGHTS
      raw.wt <- weight.binomial(x = xx, 
                                y = yy, 
                                beta = c(a00, raw.coefficients),  # Or beta_with_int
                                intercept = intercept, 
                                del = del) # del = 0.125 by standard
      
      ## STEP 5: Scaling data, now using reweighted weights (based on STEP 4.)
      sclw <- prepara(x = xx, 
                      y = yy, 
                      family = family, 
                      index = which(raw.wt == 1), # Only for those which are outlier-free
                      robu = 0)  # Indeed non-robust since this should be outlier-free as well
      xss <- sclw$xnor # Extracting standardized X
      yss <- sclw$ycen # Extracting centered y (but for Binomial no centering is actually performed)
      # NOTE: This gives us a new dataset, for which we -in essence- can tune again
      
      ## STEP 6a: FITTING ALL HYPERPARAMETERS THEIR MODELS (CV TUNING OBJECT)
      # STEP 6a: Fitting all hyperparameters their models if no lambdaw is given (full regularization path for alphas)
      if (missing(lambdaw)) {
        reweighted_cv <- cva.glmnet(x = xss[which(raw.wt == 1), ],
                                    y = yss[which(raw.wt == 1)],
                                    family = "binomial",
                                    nfolds = 10, # TODO: CHOOSE
                                    type.measure = "deviance", # This should be passable to lower lying functions such as cv.glmnet
                                    standardize = FALSE, # Since data already standardized!
                                    intercept = TRUE) # True for binomial!
        
        ## STEP 6b. CHOOSING OPTIMAL ALPHA (JUST MIN)
        alphaw_index <- which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm))) # Extracting index (bit weird but OK)
        alphaw_best <- reweighted_cv$alpha[alphaw_index] # Actual value of the best alpha!
        
        ## STEP 6c. CHOOSING OPTIMAL LAMBDA WITHIN OPTIMAL ALPHA # TODO: MIGHT DO THIS MORE INDEPENDENTLY
        if (type_lambdaw == "min") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.min # NEW!
        } else if (type_lambdaw == "1se") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.1se
        }
        
        # STEP 6a: 'Tuning' hyperparameters with if single lambdaw is given (UNLIKELY, i.e. the best lambdaw is the only one given!)
      } else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw_best <- lambdaw
        
        # STEP 6a. Tuning hyperparamters if vector of lambdaw is given (UNLIKELY)
      } else if (!missing(lambdaw) & length(lambdaw) > 1) { # Multiple lambdaw given
        reweighted_cv <- cva.glmnet(x = xss[which(raw.wt == 1), ],  # With the new reweighted dataset
                                    y = yss[which(raw.wt == 1)],  # With the new reweighted dataset
                                    family = "binomial", 
                                    lambda = lambdaw, # Passing on lambdaw vector given by user
                                    nfolds = 10, # TODO CHOOSE
                                    type.measure = "deviance",  # For binomial models negative loglikelihood
                                    standardize = FALSE, 
                                    intercept = TRUE)
        
        ## STEP 6b. CHOOSING OPTIMAL ALPHA (JUST MIN)
        alphaw_index <- which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm))) # Extracting index (bit weird but OK)
        alphaw_best <- reweighted_cv$alpha[alphaw_index] # Actual value of the best alpha!
        
        ## STEP 6c. CHOOSING OPTIMAL LAMBDA WITHIN OPTIMAL ALPHA # TODO: MIGHT DO THIS MORE INDEPENDENTLY
        if (type_lambdaw == "min") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.min # NEW!
          alphaw_best <- alphabest # NEW: Just repassin the alpha as well (even though I could retune for alpha as well)
        } else if (type_lambdaw == "1se") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.1se
        }
      }
      
      ## STEP 7: FITTING FINAL REWEIGHTED MODEL USING OPTIMAL HYPERPARAMETERS (I.E. ON REWEIGHTED DATA)
      print("alphabest:") # TODO REMOVE
      print(alphabest) # TODO REMVOE
      print("lambdaw:") # TODO REMOVE
      print(lambdaw) # TODO REMOVE
      
      # STEP 7: Binomial fitting final model using optimal hyperparameters on reweighted data
      if (family == "binomial") {
        fitw <- glmnet(x = xss[which(raw.wt == 1), ], 
                       y = yss[which(raw.wt == 1)], 
                       family = family, 
                       alpha = alphabest, 
                       lambda = lambdaw, 
                       standardize = FALSE, 
                       intercept = TRUE) # NEW: changed this to TRUE for binomial because standardization doesn't work that way for binomial
      }
      
      # TODO REMOVE
      print("I am past fitting the final reweighted fit object")
      
      ## STEP 8a. INTERCEPT HANDLING OF REWEIGHTED FIT
      if (isFALSE(intercept)) {
        a0 <- 0
      } else (isTRUE(intercept)) {
        a0 <- fitw$a0  # Just get it from the model object
      }
      
      ## STEP 8b. COEFFICIENT HANDLING OF REWEIGHTED FIT
      coefficients <- drop(as.matrix(fitw$beta)/sclw$sigx)
      
      # TODO REMOVE SOME PRINT STATEMENTS
      print("Structure of coefficients:") # TODO REMOVE
      print(str(coefficients)) # # TODO REMOVE
      print("Structure of a0:") # # TODO REMOVE
      print(str(a0)) # TODO REMOVE
      print("Calculated coefs")  # TODO REMOVE
      
      ## STEP 9 CALCULATING WEIGHTS
      wgt <- weight.binomial(x = xx, 
                             y = yy, 
                             beta = c(a0, coefficients), 
                             intercept = intercept, 
                             del = del)
      
      ## STEP 10: CALCULATING REWEIGHTED RESIDUALS (TODO WHATS DONE WITH THESE)
      reweighted.residuals <- -(yy * cbind(1, xx) %*% c(a0, coefficients)) + log(1 + exp(cbind(1, xx) %*% c(a0, coefficients)))
      
      #----------------------------- END OF MAIN (scal == TRUE) BINOMIAL CASE -----------------------#  
      #----------------------------- START OF MAIN (scal == TRUE) GAUSSIAN CASE -----------------------#  
      
      ## STEP 3. Handling of the results from the final (non-reweighted) fit.
    } else if (family == "gaussian") { # GAUSSIAN
      # STEP 3a. Intercept handling: Case Gaussian
      if (isFALSE(intercept)) {
        a00 <- 0
      } else if (isTRUE(intercept)) {
        a00 <- drop(scl$muy + fit$a0 - as.vector(as.matrix(fit$beta)) %*% (scl$mux/scl$sigx))
      }
      
      # STEP 3b. Coefficient handling: Case Gaussian
      raw.coefficients <- drop(as.matrix(fit$beta)/scl$sigx)
      
      # STEP 3c. Residuals handling
      raw.residuals <- yy - cbind(1, xx) %*% c(a00, raw.coefficients)
      raw.rmse <- sqrt(mean(raw.residuals^2))  # In-sample loss!
      
      ## STEP 4. Calculating Weights
      raw.wt <- weight.gaussian(resi = raw.residuals, 
                                ind = indexbest, 
                                del = del)$we
      
      ## STEP 5. SCALING USING THE CALCULATED WEIGHTS
      sclw <- prepara(x = xx, 
                      y = yy, 
                      family = family, 
                      index = which(raw.wt == 1), 
                      robu = 0)
      xss <- sclw$xnor  # Extracting normalized X
      yss <- sclw$ycen  # EXtracting centered y
      
      ## STEP 6. RETUNING USING THE REWEIGHTING DATASET (Gaussian)
      # STEP 6a. Fitting all models their hyperparameters (cva.glmnet()): Case: no lambdaw given (USUAL, i.e. full regularization path fitted)
      if ((missing(lambdaw))) {
        reweighted_cv <- cva.glmnet(x = xss[which(raw.wt == 1), ],
                                    y = yss[which(raw.wt == 1)], 
                                    family = "gaussian", 
                                    nfolds = 10, 
                                    standardize = FALSE, # Because xss/yss is already standardized
                                    intercept = FALSE, # FALSE OK because Gaussian, hence goes through origin 
                                    type.measure = "mse") # NEW: REMOVED $lambda.min here because we extract it later anyways
        
        # STEP 6b. Choosing optimal alpha (just min CV error alpha)
        alphaw_index <- which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm))) # Extracting index (bit weird but OK)
        alphaw_best <- reweighted_cv$alpha[alphaw_index] # Actual value of the best alpha!
        
        ## STEP 6c. CHOOSING OPTIMAL LAMBDA WITHIN OPTIMAL ALPHA # TODO: MIGHT DO THIS MORE INDEPENDENTLY (i.e. not first look at alpha)
        if (type_lambdaw == "min") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.min # NEW!
        } else if (type_lambdaw == "1se") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.1se # NEW!
        }
        
        # STEP 6a. Case: if single hyperparamters, the given one is the optimal one
      } else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw_best <- lambdaw
        alphaw_best <- alphabest # Also just assigning alpha from before # TODO a bit shortcut but will be OK
        
        # STEP 6a. Case: lambdaw vector is given: tuning using given lambdaw
      } else if (!missing(lambdaw) & length(lambdaw) > 1) {
        reweighted_cv <- cva.glmnet(x = xss[which(raw.wt == 1), ],
                                    y = yss[which(raw.wt == 1)], 
                                    family = family, 
                                    lambda = lambdaw, 
                                    nfolds = 5, 
                                    standardize = FALSE, # Because xss/yss is already standardized
                                    intercept = FALSE,  # FALSE is OK because Gaussian, hence goes through origin
                                    type.measure = "mse") # NEW: REMOVED $lambda.min here because we extract it later anyways
        
        # STEP 6b. Choosing optimal alpha (just min CV error alpha)
        alphaw_index <- which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm))) # Extracting index (bit weird but OK)
        alphaw_best <- reweighted_cv$alpha[alphaw_index] # Actual value of the best alpha!
        
        ## STEP 6c. CHOOSING OPTIMAL LAMBDA WITHIN OPTIMAL ALPHA # TODO: MIGHT DO THIS MORE INDEPENDENTLY (i.e. not first look at alpha)
        if (type_lambdaw == "min") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.min # NEW!
        } else if (type_lambdaw == "1se") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.1se # NEW!
        }
      }
      
      ## STEP 7: FITTING FINAL (REWEIGHTED) MODEL ON REWEIGHTED DATA
      fitw <- glmnet(x = xss[which(raw.wt == 1), ], 
                     y = yss[which(raw.wt == 1)], 
                     family = "gaussian", 
                     alpha = alphaw_best, 
                     lambda = lambdaw_best, 
                     standardize = FALSE, 
                     intercept = FALSE)
      
      ## STEP 8. RESULTS HANDLING OF THE REWEIGHTED FIT
      # STEP 8a. Intercept handling
      if (isFALSE(intercept)) {
        a0 <- 0
      } else if (isTRUE(intercept)) {
        a0 <- drop(sclw$muy + fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% (sclw$mux/sclw$sigx))
      }
      
      # STEP 8b. Coefficient (reweighted) handling
      coefficients <- drop(as.matrix(fitw$beta)/sclw$sigx)
      
      # STEP 8c. Residuals (reweighted) handling
      reweighted.residuals <- yy - cbind(1, xx) %*% c(a0, coefficients)
      reweighted.rmse <- sqrt(mean(reweighted.residuals^2))
      
      ## STEP 9: 
      wgt <- weight.gaussian(resi = reweighted.residuals, 
                             ind = raw.wt == 1, 
                             del = del)$we
    } # End of Case: Gaussian
    # End of isTRUE(scal)
    
    #----------------------------------SCALING FALSE --------------------------------#
  } else if (isFALSE(scal)) { # Think (!) this amounts to having isFALSE(scal):
    # Case: No scaling (unlikely)
    
    ## STEP 1. Scaling (if scal TRUE) using the best indices found (related to best hyperparameter), hence nonrobust (because outlier-free world)
    # DOES NOT HAPPEN HERE BECAUSE SCAL == FALSE
    
    ## STEP 2. Fitting of ordinary elastic net (glmnet) on all non-outlying data with best hyperparams and scaled outlier-free dataset
    fit <- glmnet(x = x[indexbest, ], 
                  y = y[indexbest, ], 
                  family = family, 
                  alpha = alphabest, 
                  lambda = lambdabest, 
                  standardize = FALSE, 
                  intercept = FALSE)
    
    ## STEP 3: RESULTS HANDLING OF THE NON-REWEIGHTED FIT
    # STEP 3a. Intercept handling
    if (family == "binomial") {
      if (isFALSE(intercept)){
        a00 <- 0
      } else if (isTRUE(intercept)){
        a00 <- fit$a0 # NEW: Since intercept does not go with that weird formula
        # a00 <- drop(fit$a0 - as.vector(as.matrix(fit$beta)) %*% (sc$mux/sc$sigx))
      }
      
      # STEP 3b. Coefficient handling
      raw.coefficients <- drop(as.matrix(fit$beta)/sc$sigx)
      
      # STEP 3c. Residuals handling
      raw.residuals <- -(y * x %*% as.matrix(fit$beta)) + log(1 + exp(x %*% as.matrix(fit$beta)))
      
      ## STEP 4. CALCULATING WEIGHTS
      raw.wt <- weight.binomial(x = xx, 
                                y = yy, 
                                beta = c(a00, raw.coefficients), 
                                intercept = intercept, 
                                del = del)
      
      ## STEP 5. Scaling (if scal TRUE) using the best REWEIGHTED indices found found (related to best hyperparamters), hence nonrobust but potentially less conservative
      # NOTE: NOTE DONE BECAUSE SCAL == FALSE
      
      ## STEP 6: HYPERPARAMETER TUNING
      # STEP 6a. Fitting all hyperparameters their models: Case: no lambdaw given (USUALLY, fits whole regularization path)
      if (missing(lambdaw)) {
        reweighted_cv <- cva.glmnet(x = x[which(raw.wt == 1), ], 
                                    y = y[which(raw.wt == 1)], 
                                    family = "binomial", 
                                    nfolds = 10,
                                    type.measure = "deviance",
                                    standardize = FALSE, 
                                    intercept = FALSE)
        
        # STEP 6b. Extracting optimal reweighted alpha (JUST TAKING ONE WHERE MINIMUM ERROR IS OBTAINED)
        alphaw_index <- which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm))) # Extracting index (bit weird but OK)
        alphaw_best <- reweighted_cv$alpha[alphaw_index] # Actual value of the best alpha!
        
        # STEP 6c. Extracting otimal reweighted lambda
        if (type_lambdaw == "min") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.min # NEW!
        } else if (type_lambdaw == "1se") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.1se # NEW!
        }
        
        # STEP 6a. Fitting all hyperparameters their models: Case: single lambdaw given, hence it is the best
      } else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw_best <- lambdaw
        alphaw_best <- alphasbest # Just passing on the best alpha from the non-reweighted
      } else if (!missing(lambdaw) & length(lambdaw) > 1) {
        reweighted_cv <- cva.glmnet(x = x[which(raw.wt == 1), ], 
                                    y = y[which(raw.wt == 1)], 
                                    family = "binomial", 
                                    nfolds = 10, 
                                    type.measure = "deviance",
                                    standardize = FALSE, 
                                    intercept = FALSE)
        
        # STEP 6b. Extracting optimal reweighted alpha (JUST TAKING ONE WHERE MINIMUM ERROR IS OBTAINED)
        alphaw_index <- which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm))) # Extracting index (bit weird but OK)
        alphaw_best <- reweighted_cv$alpha[alphaw_index] # Actual value of the best alpha!
        
        # STEP 6c. Extracting otimal reweighted lambda
        if (type_lambdaw == "min") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.min # NEW!
        } else if (type_lambdaw == "1se") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.1se # NEW!
        }
      }
      
      ## STEP 7: FINAL REWEIGHTED FITTING 
      # Fitting using optimally found alphaw_best and lambdaw_best
      fitw <- glmnet(x = x[which(raw.wt == 1), ], 
                     y = y[which(raw.wt == 1)], 
                     family = "binomial", 
                     alpha = alphaw_best, 
                     lambda = lambdaw_best, 
                     standardize = FALSE, 
                     intercept = TRUE) # BECAUSE STILL CASE IS BINOMIAL
      
      ## STEP 8: RESULTS HANDLING OF REWEIGHTED FIT
      # STEP 8a. Intercept handling of reweighted fit
      if(isFALSE(intercept)){
        a0 <- 0
      } else if (isTRUE(intercept)) {
        a0 <- fitw$a0  # NEW: CHANGING TO THIS
      }
      
      # STEP 8b. Coefficients handling of reweighted fit
      coefficients <- drop(as.matrix(fitw$beta)/sc$sigx)
      
      ## STEP 9: CALCULATING WEIGHTS AGAIN
      wgt <- weight.binomial(x = xx, 
                             y = yy, 
                             beta = c(a0, coefficients), 
                             intercept = intercept, 
                             del = del)
      
      
      reweighted.residuals <- -(yy * cbind(1, xx) %*% c(a0, coefficients)) + log(1 + exp(cbind(1, xx) %*% c(a0, coefficients)))
      #----------------------- END OF BINOMIAL WITHIN SCAL == FALSE ------------------ #
      #----------------------- BEGIN OF GAUSSIAN WITHIN SCAL == FALSE ------------------ #
    } else if (family == "gaussian") { # Case: Gaussian (and no scaling)
      
      ## STEP 1. Scaling (if scal TRUE) using the best indices found (related to best hyperparameter), hence nonrobust (because outlier-free world)
      # DOES NOT HAPPEN HERE BECAUSE SCAL == FALSE
      
      ## STEP 3. RESULTS HANDLING OF THE NON-REWEIGHTED FIT
      # STEP 3a. Intercept handling of non reweighted fit
      if (isFALSE(intercept)) {
        a00 <- 0
      }  else if (isTRUE(intercept)) {
        a00 <- drop(sc$muy + fit$a0 - as.vector(as.matrix(fit$beta)) %*% (sc$mux/sc$sigx))
      }
      
      # STEP 3b. Coefficient handling of non-reweighted fit
      raw.coefficients <- drop(as.matrix(fit$beta)/sc$sigx)
      
      # STEP 3c. Residuals handling of non-reweighted fit
      raw.residuals <- yy - cbind(1, xx) %*% c(a00, raw.coefficients)
      raw.rmse <- sqrt(mean(raw.residuals^2))
      
      ## STEP 4. CALCULATIONS OF WEIGHTS FOR REWEIGHTED DATASET
      raw.wt <- weight.gaussian(resi = raw.residuals, 
                                ind = indexbest, 
                                del = del)$we
      
      ## STEP 5. Scaling (if scal TRUE) using the best REWEIGHTED indices found found (related to best hyperparamters), hence nonrobust but potentially less conservative
      # NOTHING DONE BECAUSE scal == FALSE
      
      ## STEP 6. TUNING HYPERPARAMETER ON THE REWEIGHTED SET
      # STEP 6a. Fitting all the models their hyperparamters (i.e. call cva.glmnet())
      if (missing(lambdaw)) {
        reweighted_cv <- cva.glmnet(x = x[which(raw.wt == 1), ], 
                                    y = y[which(raw.wt == 1)], 
                                    family = "gaussian", 
                                    nfolds = 10, 
                                    standardize = FALSE, 
                                    intercept = FALSE,  # OK because Gaussian
                                    type.measure = "mse")
        # STEP 6b. Extracting optimal reweighted alpha (JUST TAKING ONE WHERE MINIMUM ERROR IS OBTAINED)
        alphaw_index <- which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm))) # Extracting index (bit weird but OK)
        alphaw_best <- reweighted_cv$alpha[alphaw_index] # Actual value of the best alpha!
        
        # STEP 6c. Extracting otimal reweighted lambda
        if (type_lambdaw == "min") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.min # NEW!
        } else if (type_lambdaw == "1se") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.1se # NEW!
        }
        # STEP 6a. If only a single lambdaw is given, automatically it is the best
      } else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw_best <- lambdaw
        alphaw_best <- alphabest # Just passing on the best alpha from before
      } else if (!missing(lambdaw) & length(lambdaw) > 1) {
        reweighted_cv <- cva.glmnet(x = x[which(raw.wt == 1), ], 
                                    y = y[which(raw.wt == 1)], 
                                    family = "gaussian", 
                                    lambda = lambdaw, 
                                    nfolds = 10, 
                                    standardize = FALSE,
                                    intercept = FALSE,  # OK because Gaussian
                                    type.measure = "mse")
        # STEP 6b. Extracting optimal reweighted alpha (JUST TAKING ONE WHERE MINIMUM ERROR IS OBTAINED)
        alphaw_index <- which.min(sapply(reweighted_cv$modlist, function(mod) min(mod$cvm))) # Extracting index (bit weird but OK)
        alphaw_best <- reweighted_cv$alpha[alphaw_index] # Actual value of the best alpha!
        
        # STEP 6c. Extracting otimal reweighted lambda
        if (type_lambdaw == "min") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.min # NEW!
        } else if (type_lambdaw == "1se") {
          lambdaw_best <- reweighted_cv$modlist[[alphaw_index]]$lambda.1se # NEW!
        }
      }
      
      ## STEP 7: FINAL REWEIGHTED FITTING 
      fitw <- glmnet(x = x[which(raw.wt == 1), ], 
                     y = y[which(raw.wt == 1)], 
                     family = "gaussian", 
                     alpha = alphabest, 
                     lambda = lambdaw, 
                     standardize = FALSE, 
                     intercept = FALSE) # DOES THIS MAKE SENSE STILLL -> NO NOT IF NOT STANDARDIZED BY ITSELF SOMEWHERE BEFORE...
      
      ## STEP 8: RESULTS HANDLING OF THE FINAL REWEIGHTED FIT (Case: scal: FALSE and family: "Gaussian")
      # STEP 8a. Intercept handling
      if (isFALSE(intercept)) {
        a0 <- 0
      } else if (isTRUE(intercept)) {
        a0 <- drop(sc$muy + fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% (sc$mux/sc$sigx))
      }
      
      # STEP 8b. Coefficient handling
      coefficients <- drop(as.matrix(fitw$beta)/sc$sigx)
      
      # STEP 8c. Residuals handling
      reweighted.residuals <- yy - cbind(1, xx) %*% c(a0, coefficients)
      reweighted.rmse <- sqrt(mean(reweighted.residuals^2))
      
      ## STEP 9. FINAL WEIGHT CALCULATION
      wgt <- weight.gaussian(resi = reweighted.residuals, 
                             ind = raw.wt == 1, 
                             del = del)$we
    } # End family = "gaussian"
  } # End scal == FALSE
  
  # Calculating number of nonzero coefficients
  num.nonzerocoef <- sum(coefficients != 0)

  ## PREPARING OUTPUT
  if(family == "binomial"){
    output <- list(indexbest = indexbest, 
                   raw.wt = raw.wt, 
                   wgt = wgt, 
                   a00 = a00, 
                   raw.coefficients = raw.coefficients, 
                   a0 = a0, 
                   coefficients = coefficients, 
                   alphabest = alphabest,
                   alphaw_best = alphaw_best,
                   lambdabest = lambdabest, 
                   lambdaw_best = lambdaw_best, 
                   raw.residuals = raw.residuals, 
                   residuals = reweighted.residuals,
                   num.nonzerocoef = num.nonzerocoef)
    
  } else if (family == "gaussian"){
    output <- list(indexbest = indexbest, 
                   raw.wt = raw.wt, 
                   wgt = wgt, 
                   a00 = a00, 
                   raw.coefficients = raw.coefficients, 
                   a0 = a0, 
                   coefficients = coefficients, 
                   alphabest = alphabest,
                   alphaw_best = alphaw_best,
                   lambdabest = lambdabest, 
                   lambdaw_best = lambdaw_best, 
                   raw.residuals = raw.residuals, 
                   residuals = reweighted.residuals, 
                   raw.rmse = raw.rmse, 
                   reweighted.rmse = reweighted.rmse,
                   num.nonzerocoef = num.nonzerocoef)
  }
}