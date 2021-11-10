interacted_matrix <- function(X, w, interacted){
  
  X.nam <- colnames(X)
  
  if(is.character(interacted)){
    logic <- interacted %in% X.nam
    idx   <- which(logic)
    stopifnot(all(logic))
  } else if(is.numeric(interacted)){
    idx <- interacted
    stopifnot(max(interacted) <= ncol(X))
  } else stop("Indices must be either numeric or character")
  
  # make the matrix
  intmat <- sapply(idx, function(j) w * X[, j, drop = FALSE] )
  out <- cbind(X, w, intmat)
  colnames(out) <- c(X.nam, "w", paste0("w.", X.nam[idx]))
  out
  
} # FOR


#' Effect model
#' 
#' @param X Matrix of fixed covariates.
#' @param status Numeric vector with a unique code for each failure type. Code 0 denotes a censored observation.
#' @param w Binary vector of treatment assignment status. Equal to 1 for treatment group and 0 for control group.
#' @param interacted Variables that are to be interacted with \code{w}. Either numeric or character.
#' @param alpha The elasticnet mixing parameter for regularization in the effect model. The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1} (default) is the lasso penalty, and \code{alpha=0} the ridge penalty.
#' @param failcode Code of status that denotes the failure type of interest. Default is one.
#' @param retained Character vector of variables that shall not be regularized. Start names with \code{w.} to indicate interacted variables.
#' @param retain_w Logical. Shall treatment assignment \code{w} be retained? Defalt is \code{TRUE}.
#' @param baseline_risk User-defined baseline risk. If \code{NULL} (default), then baseline risk if estimated with \code{(X, status)} using \code{\link{baseline_risk}}.
#' @param alpha_baseline The elasticnet mixing parameter for regularization in a potential baseline risk model to obtain baseline risk. Only applicable if \code{z = NULL}. See \code{\link{baseline_risk}} for details.
#' @param ... Additional arguments to be passed.
#' 
#' @return A \code{predmod_ordinary} object.
#' 
#' @export
effect_model <- function(X, 
                         status, 
                         w, 
                         interacted = 1:ncol(X), 
                         alpha = 1,
                         failcode = 1,
                         retained = NULL, # needs to be character!
                         retain_w = TRUE,
                         baseline_risk = NULL,
                         alpha_baseline = 1,
                         ...){
  
  # input checks
  InputChecks_NA(list(X, status, w))
  CheckInputs_X(X)
  InputChecks_W(w)
  InputChecks_equal.length2(X, status)
  
  # response needs to be binary, so recode status to be binary with 1 = failure due to cause of interest
  status_bin <- ifelse(status == failcode, 1, 0)
  
  # assign variable names if there are none
  if(is.null(colnames(X))) colnames(X) <- paste0("V", 1:ncol(X))
  
  ## get interacted variables
  if(is.null(interacted)){
    
    ## no variables are interacted
    
    # get full matrix
    X_full <- cbind(X, w = w)
    
    # prepare regressor matrix with reversed w
    X_full_rev <- cbind(X, w = ifelse(w == 1, 0, 1))
    
  } else{
    
    ## at least one variable is interacted
    
    # get full matrix
    X_full <- interacted_matrix(X = X, w = w, interacted = interacted)
    
    # prepare regressor matrix with reversed w
    X_full_rev <- interacted_matrix(X = X, w = ifelse(w == 1, 0, 1), interacted = interacted)
    
  } # IF
  
  
  # get index of W
  idx_w <- which(colnames(X_full) %in% "w")
  
  # binary penalty factor. Zero means that corresponding variable isn't shrunk
  penalty.factor <- rep(1, ncol(X_full))
  
  # optional penalization of further covariates
  if(!is.null(retained)){
    
    # error checks
    stopifnot(is.character(retained))
    stopifnot(all(retained %in% colnames(X_full)))
    
    # adjust penalty factors
    retained.idx <- which(colnames(X_full) %in% retained)
    penalty.factor[retained.idx] <- 0
    
  } # IF
  
  # optional penalization of w
  if(retain_w){
    penalty.factor[idx_w] <- 0 
  } # IF
  
  # fit models
  fits <- effect_model_fit(X = X_full, 
                           status = status_bin, 
                           alpha = alpha,
                           penalty.factor = penalty.factor)
  
  # obtain risk estimates
  risk_reg <- stats::plogis(as.numeric(
    cbind(1, X_full[,fits$kept_vars, drop = FALSE]) %*% fits$coefficients$reduced))
  
  risk_rev <- stats::plogis(as.numeric(
    cbind(1, X_full_rev[,fits$kept_vars, drop = FALSE]) %*% fits$coefficients$reduced))
  
  # calculate predicted benefits
  benefits <- get_predicted_benefits(risk_reg = risk_reg, 
                                     risk_rev = risk_rev,
                                     w = w)
  
  # fit baseline risk if necessary
  if(is.null(baseline_risk)){
    
    mod_baseline <- baseline_risk(X = X, status = status,
                                  alpha = alpha_baseline, 
                                  failcode = failcode,...)
    br <- mod_baseline$risk
    
  } else{
    br <- baseline_risk
    mod_baseline <- NULL
  } # IF
  
  # estimate concordance on baseline risk
  C_outcome_baseline <- C_outcome(y = status_bin, risk = br)
  
  # return
  return(structure(list(benefits = benefits,
                        coefficients = list(baseline = mod_baseline$coefficients,
                                            full = fits$coefficients$full,
                                            reduced = summary(fits$models$reduced)$coefficients),
                        risk = list(baseline = br,
                                    regular = risk_reg,
                                    counterfactual = risk_rev),
                        concordance = list(outcome_baseline = C_outcome_baseline,
                                           outcome = C_outcome(y = status, risk = risk_reg),
                                           benefit = C_benefit(y = status_bin, 
                                                               w = w,
                                                               pred_ben = benefits$absolute)),
                        models = list(baseline = mod_baseline, 
                                      full = fits$models$full,
                                      reduced = fits$models$reduced),
                        inputs = list(status = status, status_bin = status_bin,
                                      w = w, failcode = failcode, alpha = alpha, 
                                      alpha_baseline = alpha_baseline)
  ), 
  class = "predmod_ordinary"))
  

} # FUN


effect_model_fit <- function(X, status, alpha, penalty.factor){
  
  # fit penalized regression for full model
  model_full <- glmnet::cv.glmnet(x = X, y = status, 
                                   family = "binomial",
                                   alpha = alpha, 
                                   penalty.factor = penalty.factor)
  
  # extract retained variables
  coefs_full <- glmnet::coef.glmnet(model_full, s = "lambda.min")
  colnames(coefs_full) <- "Estimate"
  
  # get indices of kept variables 
  # the 0-th index in @i here corresponds to the intercept, so drop it
  kept.vars  <- setdiff(coefs_full@i, 0)
  
  # fit final model
  X_reduced     <- X[,kept.vars, drop = FALSE]
  model_reduced <- stats::glm(y~., data = data.frame(y = status, X_reduced),
                                   family =  stats::binomial(link = "logit"))
  
  # store coefficient matrix as sparse matrix (for consistency with glmnet)
  coefs_reduced <- Matrix::Matrix(model_reduced$coefficients, sparse = TRUE, 
                                  dimnames = list(c("(Intercept)", colnames(X_reduced)), "Estimate"))
  
  return(list(models = list(full = model_full, reduced = model_reduced),
              coefficients = list(full = coefs_full, reduced = coefs_reduced),
              kept_vars = kept.vars,
              lambda_min = model_full$lambda.min))
  
} # FUN



#' Effect model with survival
#' 
#' @param X Matrix of fixed covariates.
#' @param status Numeric vector with a unique code for each failure type. Code 0 denotes a censored observation.
#' @param time vector of failure/censoring times.
#' @param w Binary vector of treatment assignment status. Equal to 1 for treatment group and 0 for control group.
#' @param interacted Variables that are to be interacted with \code{w}. Either numeric or character.
#' @param time_eval Time at at which survival shall be evaluated.
#' @param alpha The elasticnet mixing parameter for regularization in the effect model. The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1} (default) is the lasso penalty, and \code{alpha=0} the ridge penalty.
#' @param failcode Code of status that denotes the failure type of interest. Default is one.
#' @param retained Character vector of variables that shall not be regularized. Start names with \code{w.} to indicate interacted variables.
#' @param retain_w Logical. Shall treatment assignment \code{w} be retained? Defalt is \code{TRUE}.
#' @param baseline_risk User-defined baseline risk. If \code{NULL} (default), then baseline risk if estimated with \code{(X, status, time)} using \code{\link{baseline_survival}}, as 1 - S(time_eval).
#' @param alpha_baseline The elasticnet mixing parameter for regularization in a potential baseline risk model to obtain baseline risk. Only applicable if \code{z = NULL}. See \code{\link{baseline_risk}} for details.
#' 
#' @return A \code{predmod_survival} object.
#' 
#' @export
effect_model_survival <- function(X, 
                                 status, 
                                 time,
                                 w, 
                                 interacted = 1:ncol(X), 
                                 time_eval = max(time),
                                 alpha = 1,
                                 failcode = 1,
                                 retained = NULL, # needs to be character!
                                 retain_w = TRUE,
                                 baseline_risk = NULL,
                                 alpha_baseline = 1){
  
  # input checks
  InputChecks_NA(list(X, status, w))
  CheckInputs_X(X)
  InputChecks_W(w)
  InputChecks_equal.length2(X, status)
  
  # decide whether or not competing risk modeling should be used
  cmprsk_fun <- ifelse(length(unique(status)) == 2,
                       "effect_model_fit_nocmprsk",
                       "effect_model_fit_cmprsk")
  
  # response needs to be binary, so recode status to be binary with 1 = failure due to cause of interest
  status_bin <- ifelse(status == failcode, 1, 0)
  
  # assign variable names if there are none
  if(is.null(colnames(X))) colnames(X) <- paste0("V", 1:ncol(X))
  
  ## get interacted variables
  if(is.null(interacted)){
    
    ## no variables are interacted
    
    # get full matrix
    X_full <- cbind(X, w = w)
    
    # prepare regressor matrix with reversed w
    X_full_rev <- cbind(X, w = ifelse(w == 1, 0, 1))
    
  } else{
    
    ## at least one variable is interacted
    
    # get full matrix
    X_full <- interacted_matrix(X = X, w = w, interacted = interacted)
    
    # prepare regressor matrix with reversed w
    X_full_rev <- interacted_matrix(X = X, w = ifelse(w == 1, 0, 1), interacted = interacted)
    
  } # IF
  
  
  # get index of W
  idx_w <- which(colnames(X_full) %in% "w")
  
  # binary penalty factor. Zero means that corresponding variable isn't shrunk
  penalty.factor <- rep(1, ncol(X_full))
  
  # optional penalization of further covariates
  if(!is.null(retained)){
    
    # error checks
    stopifnot(is.character(retained))
    stopifnot(all(retained %in% colnames(X_full)))
    
    # adjust penalty factors
    retained.idx <- which(colnames(X_full) %in% retained)
    penalty.factor[retained.idx] <- 0
    
  } # IF
  
  # optional penalization of w
  if(retain_w){
    penalty.factor[idx_w] <- 0 
  } # IF
  
  
  # fit models
  fits <- do.call(what = get(cmprsk_fun),
                  args = list(X = X_full, 
                              X_rev = X_full_rev, 
                              status = status_bin, 
                              time = time, 
                              alpha = alpha, 
                              penalty_factor = penalty.factor,
                              failcode = failcode))
  
  # extract failure time estimates
  fail_reg <- fits$failure$regular
  fail_rev <- fits$failure$counterfactual
  
  # calculate predicted benefits
  benefits <- get_predicted_benefits(risk_reg = fail_reg, 
                                     risk_rev = fail_rev,
                                     w = w)
  
  # get failure risk at the time of interest
  risk_reg <- 1.0 - fits$funs$regular$surv(time_eval)
  risk_rev <- 1.0 - fits$funs$counterfactual$surv(time_eval)
  
  # get benefits of risk
  benefits_risk <- get_predicted_benefits(risk_reg = risk_reg, 
                                          risk_rev = risk_rev,
                                          w = w)
  
  # fit baseline risk if necessary
  if(is.null(baseline_risk)){
    
    mod_baseline <- baseline_survival(X = X, 
                                      status = status,
                                      time = time,
                                      time_eval = time_eval, 
                                      alpha = alpha_baseline,
                                      failcode = failcode)
    br <- mod_baseline$risk
    
  } else{
    br <- baseline_risk
    mod_baseline <- NULL
  } # IF baseline
  
  # estimate concordance on 1st stage
  C_outcome_baseline <- C_outcome(y = status_bin, risk = br)
  
  # return
  return(structure(list(benefits = benefits,
                        coefficients = list(baseline = mod_baseline$coefficients,
                                            full = fits$coefficients$full,
                                            reduced = fits$coefficients$reduced),
                        risk = list(baseline = br,
                                    regular = risk_reg,
                                    counterfactual = risk_rev),
                        failure = list(regular = fail_reg,
                                       counterfactual = risk_rev),
                        concordance = list(outcome_baseline = C_outcome_baseline,
                                           outcome = C_outcome(y = status_bin, risk = risk_reg),
                                           benefit = C_benefit(y = status_bin, 
                                                               w = w,
                                                               pred_ben = benefits_risk$absolute)),
                        benefits_risk = benefits_risk,
                        models = list(baseline = mod_baseline, 
                                      full = fits$models$full,
                                      reduced = fits$models$reduced),
                        inputs = list(status = status, status_bin = status_bin,
                                      time = time, w = w, failcode = failcode, 
                                      alpha = alpha, alpha_baseline = alpha_baseline)
  ), 
  class = "predmod_survival"))
  
  
} # FUN




effect_model_fit_nocmprsk <- function(X, X_rev, status, time, alpha, 
                                      penalty_factor,
                                      failcode = 1){ # dead argument here
  
  # prepare dependent variable
  y <- survival::Surv(time = time, event = status)
  
  # fit models
  model_full <- glmnet::cv.glmnet(x = X, y = y,
                                  family = "cox",
                                  type.measure = "deviance",
                                  penalty.factor = penalty_factor,
                                  alpha = alpha) 
  
  # extract retained variables
  coefs_full <- glmnet::coef.glmnet(model_full, s = "lambda.min")
  colnames(coefs_full) <- "Estimate"
  
  # get indices of kept variables 
  # there is no intercept, so need to account for zero-indexing
  kept.vars  <- coefs_full@i + 1L
  
  # fit final model
  X_reduced     <- X[,kept.vars, drop = FALSE]
  X_reduced_rev  <- X_rev[,kept.vars, drop = FALSE]
  model_reduced <- survival::coxph(survival::Surv(time = time, event = status)~.,
                                   data = data.frame(time, status, X_reduced))
  
  # obtain risk estimates
  # get linear predictor with regular and reversed treatment assignment 
  lp_reg  <- as.numeric(X_reduced %*% model_reduced$coefficients)
  lp_rev  <- as.numeric(X_reduced_rev %*% model_reduced$coefficients)
  
  # get the survival functions
  surv_obj_reg <- survival(time = time, status = status, lp = lp_reg, center = FALSE)
  surv_obj_rev <- survival(time = time, status = status, lp = lp_rev, center = FALSE)
  
  # get sorted unique failure times to obtain survival curves
  time_unique    <- sort(unique(time), decreasing = FALSE)
  surv_curve_reg <- sapply(time_unique, function(t) surv_obj_reg$surv(t))
  surv_curve_rev <- sapply(time_unique, function(t) surv_obj_rev$surv(t))
  
  # use survival curves to estimate potential failure times
  fail_reg <- expected_survival(S.hat = surv_curve_reg, Y.grid = time_unique)
  fail_rev <- expected_survival(S.hat = surv_curve_rev, Y.grid = time_unique)
  
  # coefficient summary
  temp <- summary(model_reduced)$coefficients
  est  <- temp[, "coef"]
  se   <- temp[, "se(coef)"]
  t    <- est / se
  p    <-  2 * stats::pnorm(abs(t), lower.tail = FALSE)
  cf   <- cbind(est, se, t, p)
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  # return
  return(list(models = list(full = model_full, reduced = model_reduced),
              coefficients = list(full = coefs_full, 
                                  reduced = cf),
              failure = list(regular = fail_reg, counterfactual = fail_rev),
              funs = list(regular = surv_obj_reg, counterfactual = surv_obj_rev),
              lambda_min = model_full$lambda.min))
  
} # FUN


effect_model_fit_cmprsk <- function(X, X_rev, 
                                    status, 
                                    time, alpha,
                                    penalty_factor,
                                    failcode = 1){
  
  # get the lambda grid TODO: make user choose path
  lambda_path <- get_lambda_path(x = X, time = time, 
                                 status = status, 
                                 failcode = failcode, 
                                 alpha = alpha,
                                 m = 100)
  
  # fit competing risk model with elastic net penalty along the lambda path
  model_full <- fastcmprsk::fastCrrp(fastcmprsk::Crisk(time, status, failcode = failcode) ~.,
                                     data = data.frame(time, status, X),
                                     penalty = "ENET",
                                     alpha = alpha,
                                     lambda = lambda_path,
                                     penalty.factor = penalty_factor)
  
  # get unexported function 'AIC.fcrrp'
  fnc <- utils::getFromNamespace("AIC.fcrrp", "fastcmprsk")
  
  # choose the penalty that minimizes AIC
  s <- unname(which.min(fnc(model_full)))
  
  # store it
  lambda_min <- lambda_path[s]
  
  # store coefficient matrix as sparse matrix (for consistency with glmnet)
  coefs_full <- Matrix::Matrix(model_full$coef[,s,drop=FALSE], sparse = TRUE, 
                               dimnames = list(colnames(X), "Estimate"))
  
  # get indices of kept variables 
  # there is no intercept, so need to account for zero-indexing
  kept.vars  <- coefs_full@i + 1L
  
  # fit final model
  X_reduced     <- X[,kept.vars, drop = FALSE]
  X_reduced_rev <- X_rev[,kept.vars, drop = FALSE]
  model_reduced <- cmprsk::crr(ftime = time, 
                               fstatus = status, 
                               cov1 = X_reduced, 
                               variance = TRUE,
                               failcode = failcode)
  
  # obtain risk estimates
  # get linear predictor with regular and reversed treatment assignment 
  lp_reg  <- as.numeric(X_reduced %*% model_reduced$coef)
  lp_rev  <- as.numeric(X_reduced_rev %*% model_reduced$coef)
  
  # prepare the predict object
  prep_pred <- prep_predict(time = time, status = status, k = failcode)
  
  # get the survival functions
  surv_obj_reg <- survival_cmprsk(time = time, status = status, 
                                  lp = lp_reg, prep_predict_object = prep_pred, 
                                  failcode = failcode)
  surv_obj_rev <- survival_cmprsk(time = time, status = status, 
                                  lp = lp_rev, prep_predict_object = prep_pred, 
                                  failcode = failcode)
  
  # get sorted unique failure times to obtain survival curves
  time_unique    <- sort(unique(time), decreasing = FALSE)
  surv_curve_reg <- sapply(time_unique, function(t) surv_obj_reg$surv(t))
  surv_curve_rev <- sapply(time_unique, function(t) surv_obj_rev$surv(t))
  
  # use survival curves to estimate potential failure times
  fail_reg <- expected_survival(S.hat = surv_curve_reg, Y.grid = time_unique)
  fail_rev <- expected_survival(S.hat = surv_curve_rev, Y.grid = time_unique)
  
  # coefficient summary
  est <- model_reduced$coef
  se  <- sqrt(diag(model_reduced$var))
  t   <- est / se
  p   <-  2 * stats::pnorm(abs(t), lower.tail = FALSE)
  cf  <- cbind(est, se, t, p)
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  # return
  return(list(models = list(full = model_full, reduced = model_reduced),
              coefficients = list(full = coefs_full, 
                                  reduced = cf),
              failure = list(regular = fail_reg, counterfactual = fail_rev),
              funs = list(regular = surv_obj_reg, counterfactual = surv_obj_rev),
              lambda_min = lambda_min))
  
} # FUN