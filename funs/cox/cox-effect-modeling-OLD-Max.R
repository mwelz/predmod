source(paste0(getwd(), "/funs/c-statistics/c-statistics.R"))
source(paste0(getwd(), "/funs/linear-models/effect-modeling.R"))
source(paste0(getwd(), "/funs/cox/cox-risk-modeling.R"))


#' compute log likelihood of a Cox PH model
#' 
#' @param beta vector of coefficients (p-dimensional)
#' @param time vector of time (n-dimensional)
#' @param status binary vector of status (n-dimensional)
#' @param X design matrix (nxp-dimensional)
#' 
#' @export
loglik.coxph <- function(beta, time, status, X){
  
  X     <- as.matrix(X)
  xb    <- as.numeric(X %*% beta)
  theta <- exp(xb)
  
  # logL: sum_{i:status=1} <x_i, beta> - ln(sum_{j:t_j >= t_i} theta_j )
  sum(sapply(which(status == 1), function(i ) xb[i] - log(sum(theta[time >= time[i]]))))
  
} # FUN


#' Performs variable selection based on post-selection hypothesis tests. Function follows the strategy of Wasserman and Roeder (2009, Annals of Statistics):
#' 
#' step 0.1: create matrix Z = (w, X, interaction.terms). 
#' 
#' step 0.2: partition the observations in three roughly equally-sized sets sets: D1, D2, D3.
#' 
#' step 1.1: use D1 to perform variable selection by penalized logistic regression for given value of lambda. Let Z_lambda be the design matrix associated with the retained variables.
#' 
#' step 1.2: use the retained variables from 1.1 to calculate the logistic least squares estimator beta_lambda on D1.
#' 
#' step 1.3 calculate the empirical cross-entropy loss on D2, by using beta_lambda
#' 
#' step 2: repeat steps 1.1 to 1.3 for many lambdas. Let beta_lambdahat be the logistic least squares estimator that corresponds to the retained variables of the model with the lambda that minimizes the loss. Calculate this least squares estimator on D3.
#' 
#'  step 3: perform standard hypothesis tests to decide which of the retained variables in step 2 make it to the final model. Note that the critical value needs to be adjusted; see  Wasserman and Roeder (2009, Annals of Statistics) for details
#' 
#' @param X design matrix, can also be a data frame
#' @param y vector of binary responses. 
#' @param w vector of binary treatment assignments
#' @param interacted.variables string array of variables in _X_ that shall be interacted with treatment _w_
#' @param alpha the alpha as in 'glmnet'. Default is 1, which corresponds to the Lasso
#' @param retained.variables string array of variable names in Z that shall always be retained (i.e. also interaction terms can be considered). If NULL (default), then no restriction applies. Note that treatment assignment w will always be retained by the function.
#' @param significance.level for the hypothesis tests. Default is 0.05
#' @param prediction.timeframe
#' @param lifeyears
#' 
#' @export
variable.selection.cox <- function(X, y, w,
                                   interacted.variables,
                                   lifeyears,
                                   alpha = 1,
                                   retained.variables = NULL,
                                   significance.level = 0.05){
  
  ## prepare the design matrix Z 
  n    <- nrow(X)
  p    <- ncol(X)
  
  if(is.null(colnames(X))) colnames(X) <- paste0("X", 1:p)
  
  # get design matrix for the modeling
  if(is.null(interacted.variables)){
    Z <- data.frame(w = w, X)
  } else{
    Z <- data.frame(w = w, X,
                    get.interaction.terms.matrix(X = X, 
                                                 w = w, 
                                                 interacted.variables = interacted.variables))
  } # IF
  
  # must be a matrix
  Z <- as.matrix(Z)
  
  # get Z index of forcefully retained variables
  retained.variables.Z.idx <- get.Z.index.of.retained.variables(retained.variables = retained.variables, 
                                                                X = X, Z = Z)
  
  # randomly split data into three roughly equally sized samples
  set1 <- sample(1:n, floor(n / 3), replace = FALSE)
  set2 <- sample(setdiff(1:n, set1), floor(n / 3), replace = FALSE)
  set3 <- setdiff(1:n, c(set1, set2))
  
  # prepare cross-validation: loop over the lambdas of the glmnet path
  lambdapath <- get.lambdapath(x = Z[set1,], y = y[set1])
  
  suite <- lapply(1:length(lambdapath), function(...) list(retained.variables = NA, lasso.estimates = NA))
  loss  <- rep(NA_real_, length(lambdapath))
  
  # get survival matrix; to be passed as response of glmnet
  survival.matrix_set1 <- survival::Surv(lifeyears[set1], y[set1])
  survival.matrix_set2 <- survival::Surv(lifeyears[set2], y[set2])
  survival.matrix_set3 <- survival::Surv(lifeyears[set3], y[set3])
  colnames(survival.matrix_set1) <- 
    colnames(survival.matrix_set2) <- colnames(survival.matrix_set3) <- c("time", "status")
  
  for(i in 1:length(lambdapath)){
    
    # penalized logistic regression on first set
    mod <- glmnet::glmnet(x = Z[set1,], y = survival.matrix_set1, 
                          family = "cox", type.measure = "C", 
                          alpha = alpha,
                          lambda = lambdapath[i])
    
    # obtain retained variables. Cox doesn't have intercept, so account for zero indexing
    kept.vars     <- glmnet::coef.glmnet(mod)@i + 1
    
    # make sure that all variables that are forced to be retained will be retained
    kept.vars <- unique(sort(c(kept.vars, retained.variables.Z.idx), decreasing = FALSE))
    
    # add to suite
    suite[[i]]$retained.variables <- kept.vars
    suite[[i]]$lasso.estimates    <- as.numeric(mod$beta) # Cox doesn't have intercept
    
    #  fit Cox model with retained variables on first set
    mod <- survival::coxph(survival.matrix_set1 ~., data = data.frame(Z[set1, kept.vars]))
    
    
    # calculate loss of the model on 2nd set (loss = negative logL)
    loss[i] <- -loglik.coxph(beta = mod$coefficients, 
                             time = survival.matrix_set2[,"time"], 
                             status = survival.matrix_set2[,"status"], 
                             X = Z[set2, kept.vars])
    
  } # FOR
  
  # find the lambda that minimizes the empirical loss and its associated model (S.hat)
  lambda.min   <- lambdapath[which.min(loss)]
  S.hat        <- suite[[which.min(loss)]]$retained.variables
  Z.lambda.min <- Z[, S.hat]
  
  # on third set: use S.hat to calculate Cox model with Z.lambda.min
  mod <- survival::coxph(survival.matrix_set3 ~., data = data.frame(Z.lambda.min[set3,]))
  
  # on third set: hypothesis testing
  coeffs <- summary(mod)$coefficients
  coeffs <- coeffs[, c("coef", "se(coef)", "z", "Pr(>|z|)")]
  coeffs <- cbind(coeffs, rep(NA_real_, length(S.hat)))
  colnames(coeffs) <- c("Estimate", "Std. Error", "z value", "critical value", "retain?")
  critval <- qnorm(significance.level / 2 * length(S.hat), lower.tail = FALSE)
  coeffs[, "critical value"] <- critval
  coeffs[, "retain?"] <- 1 * (abs(coeffs[, "z value"]) > critval)
  
  # final model
  D.hat   <- which(coeffs[-1, "retain?"] == 1)
  
  # make sure that all variables that are forced to be retained will be retained
  D.hat <- unique(sort(c(D.hat, retained.variables.Z.idx), decreasing = FALSE))
  
  # get the regularized estimates at the minimizing lambda
  regularized.estimates_lambda.min <- suite[[which.min(loss)]]$lasso.estimates
  names(regularized.estimates_lambda.min) <- colnames(Z)
  
  # organize output
  model.selection <- list(selected.variables = colnames(Z)[D.hat],
                          design.matrix_selected.variables = Z[, D.hat],
                          formula.final.model = paste0("Surv(lifeyears, y) ~ ", 
                                                       paste(colnames(Z)[D.hat], collapse = " + ")),
                          final.model.object = mod,
                          lambda.min = lambda.min,
                          significance.tests_lambda.min = coeffs,
                          regularized.estimates_lambda.min = regularized.estimates_lambda.min,
                          S.hat = S.hat)
  partitioning.membership <- list(set1 = sort(set1, decreasing = FALSE),
                                  set2 = sort(set2, decreasing = FALSE),
                                  set3 = sort(set3, decreasing = FALSE))
  
  return(list(final.model = model.selection, 
              partitioning.membership = partitioning.membership))
} # FUN


#' Performs effect modeling. Final effect model is done via variable selection based on post-selection hypothesis tests. Function follows the strategy of Wasserman and Roeder (2009, Annals of Statistics):
#' 
#' step 0.1: create matrix Z = (w, X, interaction.terms). 
#' 
#' step 0.2: partition the observations in three roughly equally-sized sets sets: D1, D2, D3.
#' 
#' step 1.1: use D1 to perform variable selection by penalized logistic regression for given value of lambda. Let Z_lambda be the design matrix associated with the retained variables.
#' 
#' step 1.2: use the retained variables from 1.1 to calculate the logistic least squares estimator beta_lambda on D1.
#' 
#' step 1.3 calculate the empirical cross-entropy loss on D2, by using beta_lambda
#' 
#' step 2: repeat steps 1.1 to 1.3 for many lambdas. Let beta_lambdahat be the logistic least squares estimator that corresponds to the retained variables of the model with the lambda that minimizes the loss. Calculate this least squares estimator on D3.
#' 
#'  step 3: perform standard hypothesis tests to decide which of the retained variables in step 2 make it to the final model. Note that the critical value needs to be adjusted; see  Wasserman and Roeder (2009, Annals of Statistics) for details
#' 
#' @param X design matrix, can also be a data frame
#' @param y vector of binary responses. 
#' @param w vector of binary treatment assignments
#' @param interacted.variables string array of variables in _X_ that shall be interacted with treatment _w_
#' @param alpha the alpha as in 'glmnet'. Default is 1, which corresponds to the Lasso
#' @param lifeyears vector of life years. Default is NULL.
#' @param prediction.timeframe vector of the prediction time frame. Default is NULL.
#' @param retained.variables string array of variable names in Z that shall always be retained (i.e. also interaction terms can be considered). If NULL (default), then no restriction applies. Note that treatment assignment w will always be retained by the function.
#' @param significance.level for the hypothesis tests. Default is 0.05
#' 
#' @export
cox.effect.modeling <- function(X, y, w,
                                interacted.variables,
                                alpha = 1,
                                lifeyears, 
                                prediction.timeframe,
                                retained.variables = NULL,
                                significance.level = 0.05){
  
  # truncate y if necessary
  y.orig    <- y
  lifeyears <- ifelse(lifeyears <= prediction.timeframe, lifeyears, prediction.timeframe) 
  y         <- ifelse(lifeyears <= prediction.timeframe, y, 0)
  
  # make X a matrix
  X <- as.matrix(X)
  
  ### 1. fit baseline risk via Cox modeling  ----
  # no information on w allowed, so we cannot use the retained variables from the effect modeling
  baseline.mod  <- baseline.risk.cox(X = X, y = y, alpha = alpha, 
                                     lifeyears = lifeyears,
                                     prediction.timeframe = prediction.timeframe)
  baseline.risk <- baseline.mod$response 
  
  ### 2. perform variable selection ----
  # We use the strategy in Wasserman and Roeder (2009; Annals of Statistics)
  vs.obj <- variable.selection.cox(X = X, y = y, w = w, 
                                   lifeyears = lifeyears,
                                   interacted.variables = interacted.variables, 
                                   alpha = alpha, 
                                   retained.variables = retained.variables, 
                                   significance.level = significance.level)
  
  # get design matrix associated with the final model
  Z     <- vs.obj$final.model$design.matrix_selected.variables
  
  ### 3. fit the final effect model and use it for risk estimates ----
  # fit the final model, this time on full sample (no penalty required, selection already took place!)
  final.model        <- survival::coxph(survival::Surv(lifeyears, y) ~., data = as.data.frame(Z)) 
  predicted.benefits <- effect.model.predicted.benefits(X = X, y = y, w = w, 
                                                        final.model = final.model, Z = Z)
  
  # coefficients
  coefficients <- summary(final.model)$coefficients
  
  ### 4. return ----
  return(list(
    inputs = list(X = X, w = w, y = y.orig, 
                  lifeyears = lifeyears, 
                  prediction.timeframe = prediction.timeframe, 
                  y.prediction.timeframe = y),
    average.treatment.effect = mean(predicted.benefits$pred.ben.abs),
    baseline.model = baseline.mod,
    effect.model = list(formula = vs.obj$final.model$formula.final.model, 
                        selected.data = Z,
                        glm.obj = final.model,
                        coefficients = coefficients[,1],
                        summary = coefficients,
                        model.building = vs.obj),
    risk = list(risk.regular.w = predicted.benefits$risk.regular.w,
                risk.flipped.w = predicted.benefits$risk.flipped.w,
                risk.baseline  = baseline.risk),
    benefits = list(predicted.absolute.benefit = predicted.benefits$pred.ben.abs,
                    predicted.relative.benefit = predicted.benefits$pred.ben.rel,
                    predicted.absolute.benefit.raw = predicted.benefits$pred.ben.abs.raw,
                    predicted.relative.benefit.raw = predicted.benefits$pred.ben.rel.raw),
    C.statistics = list(c.index.outcome = C.index.outcome(y = y, 
                                                          risk.prediction = predicted.benefits$risk.regular.w),
                        c.index.benefit = C.index.benefit(y = y, w = w, 
                                                          predicted.benefit = predicted.benefits$pred.ben.abs))
  ))
  
} # FUN
