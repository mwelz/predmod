#' Performs risk modeling with a competing risk model.
#' 
#' The linear predictor of the competing risk model is defined as 
#' \deqn{\alpha + \beta_1 w + \beta_2 z + \beta_3 z * w ,}
#' where \eqn{z} is the typically equal to the linear predictor of a baseline risk model. If \code{constant = TRUE}, it is enforced that \eqn{\beta_3 = 0}. 
#' 
#' @param X Matrix of fixed covariates.
#' @param status Numeric vector with a unique code for each failure type. Code 0 denotes a censored observation.
#' @param time vector of failure/censoring times.
#' @param w Binary vector of treatment assignment status. Equal to 1 for treatment group and 0 for control group.
#' @param z The \code{z} to be used in the logistic regression model above. If \code{NULL} (default), then the linear predictor of a baseline risk model is used as \code{z}.
#' @param time_eval Time at at which survival shall be evaluated.
#' @param alpha The elasticnet mixing parameter for regularization in a potential baseline risk model to obtain the linear predictor to be used as \code{z}. Only applicable if \code{z = NULL}. See \code{\link{baseline_risk}} for details.
#' @param failcode Code of status that denotes the failure type of interest. Default is one.
#' @param constant Shall \eqn{\beta_3 = 0} be enforced in the logistic regression model above? If \code{TRUE}, then this is equivalent to assuming that there is a constant treatment effect. Default is \code{FALSE}.
#' @param ... Additional arguments to be passed.
#' 
#' @return A \code{predmod_ordinary} object.
#' 
#' @export
risk_model_survival <- function(X, 
                                status, 
                                time,
                                w, 
                                z = NULL, 
                                time_eval = max(time),
                                alpha = 1,
                                failcode = 1,
                                constant = FALSE,
                                ...)
{
  
  # input checks
  InputChecks_NA(list(X, status))
  CheckInputs_X(X)
  InputChecks_equal.length3(X, status, time)
  stopifnot(length(time_eval) == 1)
  
  # assign variable names if there are none
  if(is.null(colnames(X))) colnames(X) <- paste0("V", 1:ncol(X))
  
  # decide whether or not competing risk modeling should be used
  cmpr <- ifelse(length(unique(status)) == 2, 
                 yes = "risk_model_survival_nocmprisk", 
                 no = "risk_model_survival_cmprisk")
  
  # call correct main function
  do.call(what = get(cmpr),
          args = list(X = X, 
                      status = status, 
                      time = time,
                      w = w, 
                      z = z, 
                      time_eval = time_eval,
                      alpha = alpha,
                      failcode = failcode,
                      constant = constant,
                      ... = ...))
  
} # FUN



risk_model_survival_nocmprisk <- function(X, 
                                          status, 
                                          time,
                                          w, 
                                          z = NULL, 
                                          time_eval = max(time),
                                          alpha = 1,
                                          failcode = 1,
                                          constant = FALSE,
                                          ...){
  
  # input checks
  InputChecks_NA(list(X, status, w))
  CheckInputs_X(X)
  InputChecks_W(w)
  InputChecks_equal.length2(X, status)
  
  # assign variable names if there are none
  if(is.null(colnames(X))) colnames(X) <- paste0("V", 1:ncol(X))
  
  # make response binary for later use in concordance, with 1 = failure due to cause of interest
  status_bin <- ifelse(status == failcode, 1, 0)
  
  ## stage 1
  if(is.null(z)){
    
    # in case no z is supplied, fit a baseline risk model
    stage1 <- baseline_survival_nocmprsk(X = X, 
                                         status = status,
                                         time = time,
                                         time_eval = time_eval, 
                                         alpha = alpha,
                                         failcode = failcode, 
                                         ... = ...)
    
    # take the baseline model's linear predictor as z
    z <- stage1$linear_predictor
    
    # extract baseline risk
    baseline.risk <- stage1$risk
    
  } else{
    
    # in case values for z are supplied, don't fit baseline risk model
    stage1  <- baseline.risk <- NULL
    
  } # IF
  
  
  ## stage 2
  stage2 <- risk_model_stage2_nocmprsk(status = status_bin, 
                                       time = time, 
                                       w = w, z = z, 
                                       w_flipped = ifelse(w == 1, 0, 1), 
                                       constant = constant)
  
  # get failure risk at the time of interest
  risk_reg <- 1.0 - stage2$funs$regular$surv(time_eval)
  risk_rev <- 1.0 - stage2$funs$counterfactual$surv(time_eval)
  
  # get benefits of risk
  benefits_risk <- get_predicted_benefits(risk_reg = risk_reg, 
                                          risk_rev = risk_rev,
                                          w = w)
  
  # coefficient summary
  temp <- get_coefs(stage2$model)
  est  <- temp[, "coef"]
  se   <- temp[, "se(coef)"]
  t    <- est / se
  p    <-  2 * stats::pnorm(abs(t), lower.tail = FALSE)
  cf   <- cbind(est, se, t, p)
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  

  # return
  return(structure(list(benefits = benefits_risk,
                        coefficients = list(baseline = stage1$coefficients,
                                            stage2 = cf),
                        risk = list(baseline = baseline.risk,
                                    regular = risk_reg,
                                    counterfactual = risk_rev),
                        funs = stage2$funs,
                        models = list(baseline = stage1$model, stage2 = stage2$model),
                        inputs = list(status = status, status_bin = status_bin,
                                      time = time, time_eval = time_eval, 
                                      w = w, failcode = failcode, z = z, 
                                      constant = constant, alpha = alpha)
  ), 
  class = "risk_model_surv"))
  
  
} # FUN



risk_model_stage2_nocmprsk <- function(status, time, w, z,
                                        w_flipped = ifelse(w == 1, 0, 1), 
                                        constant = FALSE,
                                        ...)
{
  if(constant){
    
    # prepare X for second stage...
    X_stage2 <- cbind(w = w, z = z)
    
    # ... and its counterpart with flipped w
    X_stage2_rev <- cbind(w = w_flipped, z = z)
    
  } else{
    
    # prepare X for second stage...
    X_stage2 <- cbind(w = w, z = z, w.z = w * z)
    
    # ... and its counterpart with flipped w
    X_stage2_rev <- cbind(w = w_flipped, z = z, w.z = w_flipped * z)
    
  } # IF
  
  
  # fit the model
  model <- survival::coxph(survival::Surv(time = time, event = status)~.,
                           data = data.frame(time, status, X_stage2), ...)
  
  # get the coefficients and account for possibility of NA-coefs (due to multicollinearity)
  cfs  <- model$coefficients
  keep <- !is.na(cfs)
  
  # get linear predictor with regular and reversed treatment assignment 
  lp_reg  <- as.numeric(X_stage2[,keep, drop = FALSE] %*% model$coefficients[keep])
  lp_rev  <- as.numeric(X_stage2_rev[,keep, drop = FALSE] %*% model$coefficients[keep])
  
  # get the survival functions
  surv_obj_reg <- survival(time = time, status = status, lp = lp_reg, center = FALSE)
  surv_obj_rev <- survival(time = time, status = status, lp = lp_rev, center = FALSE)
  
  # return
  return(list(model = model,
              funs = list(regular = surv_obj_reg, counterfactual = surv_obj_rev)
  ))
} # FUN



risk_model_survival_cmprisk <- function(X, 
                                        status, 
                                        time,
                                        w, 
                                        z = NULL, 
                                        time_eval = max(time),
                                        alpha = 1,
                                        failcode = 1,
                                        constant = FALSE,
                                        ...){
  
  # input checks
  InputChecks_NA(list(X, status, w))
  CheckInputs_X(X)
  InputChecks_W(w)
  InputChecks_equal.length2(X, status)
  
  # assign variable names if there are none
  if(is.null(colnames(X))) colnames(X) <- paste0("V", 1:ncol(X))
  
  # make response binary for later use in concordance, with 1 = failure due to cause of interest
  status_bin <- ifelse(status == failcode, 1, 0)
  
  ## stage 1
  if(is.null(z)){
    
    # in case no z is supplied, fit a baseline risk model
    stage1 <- baseline_survival_cmprsk(X = X, 
                                       status = status,
                                       time = time,
                                       time_eval = time_eval,
                                       alpha = alpha,
                                       failcode = failcode,
                                       ... = ...)
    
    # take the baseline model's linear predictor as z
    z <- stage1$linear_predictor
    
    # extract baseline risk
    baseline.risk <- stage1$risk
    
  } else{
    
    # in case values for z are supplied, don't fit baseline risk model
    stage1 <- baseline.risk <- NULL
    
  } # IF
  
  
  ## stage 2
  stage2 <- risk_model_stage2_cmprsk(status = status_bin, 
                                     time = time, 
                                     w = w, z = z, 
                                     w_flipped = ifelse(w == 1, 0, 1), 
                                     constant = constant,
                                     failcode = failcode)

  # get failure risk at the time of interest
  risk_reg <- 1.0 - stage2$funs$regular$surv(time_eval)
  risk_rev <- 1.0 - stage2$funs$counterfactual$surv(time_eval)
  
  # get benefits of risk
  benefits_risk <- get_predicted_benefits(risk_reg = risk_reg, 
                                          risk_rev = risk_rev,
                                          w = w)
  
  # coefficient summary
  est <- stage2$model$coef
  se  <- sqrt(diag(stage2$model$var))
  t   <- est / se
  p   <-  2 * stats::pnorm(abs(t), lower.tail = FALSE)
  cf  <- cbind(est, se, t, p)
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  
  # return
  return(structure(list(benefits = benefits_risk,
                        coefficients = list(baseline = stage1$coefficients,
                                            stage2 = cf),
                        risk = list(baseline = baseline.risk,
                                    regular = risk_reg,
                                    counterfactual = risk_rev),
                        funs = stage2$funs,
                        models = list(baseline = stage1$model, stage2 = stage2$model),
                        inputs = list(status = status, status_bin = status_bin,
                                      time = time, time_eval = time_eval,
                                      w = w, failcode = failcode, z = z, 
                                      constant = constant, alpha = alpha)
  ), 
  class = "risk_model_surv"))
  
} # FUN



risk_model_stage2_cmprsk <- function(status, time, w, z,
                                     w_flipped = ifelse(w == 1, 0, 1), 
                                     constant = FALSE,
                                     failcode = 1,
                                     ...)
{
  if(constant){
    
    # prepare X for second stage...
    X_stage2 <- cbind(w = w, z = z)
    
    # ... and its counterpart with flipped w
    X_stage2_rev <- cbind(w = w_flipped, z = z)
    
  } else{
    
    # prepare X for second stage...
    X_stage2 <- cbind(w = w, z = z, w.z = w * z)
    
    # ... and its counterpart with flipped w
    X_stage2_rev <- cbind(w = w_flipped, z = z, w.z = w_flipped * z)
    
  } # IF
  
  # fit the model
  model <- cmprsk::crr(ftime = time, 
                       fstatus = status, 
                       cov1 = X_stage2, variance = TRUE,
                       failcode = failcode,... = ...)
  
  # get linear predictor with regular and reversed treatment assignment 
  lp_reg  <- as.numeric(X_stage2 %*% model$coef)
  lp_rev  <- as.numeric(X_stage2_rev %*% model$coef)
  
  # prepare the predict object
  prep.pred <- prep_predict(time = time, status = status, k = failcode)
  
  # get the survival functions
  surv_obj_reg <- survival_cmprsk(time = time, status = status, 
                                  lp = lp_reg, prep_predict_object = prep.pred, 
                                  failcode = failcode)
  surv_obj_rev <- survival_cmprsk(time = time, status = status, 
                                  lp = lp_rev, prep_predict_object = prep.pred, 
                                  failcode = failcode)
  
  # return
  return(list(model = model,
              funs = list(regular = surv_obj_reg, counterfactual = surv_obj_rev)
  ))
} # FUN
