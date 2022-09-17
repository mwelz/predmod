#' baseline risk prediction via penalized logistic regression (treatment assignment is not used)
#' 
#' @param X Matrix of fixed covariates.
#' @param status Numeric vector with a unique code for each failure type. Code 0 denotes a censored observation.
#' @param alpha The elasticnet mixing parameter for regularization, with \eqn{0\le\alpha\le 1}.
#' The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1} (default) is the lasso penalty, and \code{alpha=0} the ridge penalty. If \code{NULL}, no regularization is used.
#' @param failcode Code of status that denotes the failure type of interest. Default is one.
#' @param glm_data Shall data in glm object be returned?
#' @param ... Additional arguments to be passed
#' 
#' @export
baseline_risk <- function(X,
                          status, 
                          alpha = 1,
                          failcode = 1,
                          glm_data = FALSE,
                          ...)
{
  
  # input checks
  InputChecks_NA(list(X, status))
  CheckInputs_X(X)
  InputChecks_equal.length2(X, status)
  
  # response needs to be binary, so recode status to be binary with 1 = failure due to cause of interest
  status_bin <- ifelse(status == failcode, 1, 0)
  
  # assign variable names if there are none
  if(is.null(colnames(X))) colnames(X) <- paste0("V", 1:ncol(X)) 
  
  if(is.null(alpha)){
    
    ### case 1: no regularization
    # there is no full model
    model_full <- coefs_full <- NULL
    kept_vars  <- seq_len(ncol(X))
    
  } else{
    
    ### case 2: regularization
    
    ## full model
    # fit model
    model_full <- glmnet::cv.glmnet(x = X, y = status_bin, 
                                    family = "binomial", 
                                    alpha = alpha, ...)
    
    # get coefficients at best lambda
    coefs_full <- glmnet::coef.glmnet(model_full, s = "lambda.min")
    colnames(coefs_full) <- "Estimate"
    
    # get indices of kept variables (zero-th index is intercept)
    kept_vars <- coefs_full@i 
    
    if(0L %in% kept_vars){
      kept_vars <- kept_vars[-which(kept_vars == 0L)]
    } # IF

  } # IF
  
  ## reduced model
  # fit the model
  model_red <- stats::glm(status_bin~., 
                          family = stats::binomial(link = "logit"), 
                          data = data.frame(status_bin, 
                                            X[,kept_vars,drop = FALSE]), 
                          model = FALSE, x = FALSE, y = FALSE, ...)
  
  # get coefficients and linear predictor
  coefs_red <- stats::coefficients(summary(model_red))
  lp        <- as.numeric(model_red$linear.predictors)
  
  if(!glm_data) model_red$data <- NULL
  
  # return
  return(structure(list(
    risk = as.matrix(stats::plogis(lp)),
    linear_predictor = lp,
    coefficients = list(full = coefs_full, reduced = coefs_red),
    model = list(full = model_full, reduced = model_red),
    inputs = list(status = status, 
                  status_bin = status_bin,
                  failcode = failcode,
                  alpha = alpha, 
                  X = X)
  ), class = "baseline_crss"))
  
} # FUN



#' baseline risk prediction via penalized logistic regression (treatment assignment is not used)
#' 
#' @param X Matrix of fixed covariates.
#' @param status Numeric vector with a unique code for each failure type. Code 0 denotes a censored observation.
#' @param time vector of failure/censoring times.
#' @param time_eval Time at at which survival shall be evaluated.
#' @param alpha The elasticnet mixing parameter for regularization, with \eqn{0\le\alpha\le 1}.
#' The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1} (default) is the lasso penalty, and \code{alpha=0} the ridge penalty. If \code{NULL}, no regularization is used.
#' @param failcode Code of status that denotes the failure type of interest.
#' @param ... Additional arguments to be passed
#' 
#' @export
baseline_survival <- function(X, 
                              status,
                              time,
                              time_eval = max(time),
                              alpha = 1, 
                              failcode = 1, ...)
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
                 yes = "baseline_survival_nocmprsk", 
                 no = "baseline_survival_cmprsk")
  
  # call correct main function
  do.call(what = get(cmpr),
          args = list(X = X, 
                      status = status,
                      time = time,
                      time_eval = time_eval,
                      alpha = alpha,
                      failcode = failcode,
                      ... = ...))
  
} # FUN


baseline_survival_cmprsk <- function(X, 
                                     status,
                                     time,
                                     time_eval = max(time),
                                     alpha = 1,
                                     failcode = 1,
                                     ...){
  
  if(!is.null(alpha)){
    
    # get the lambda grid TODO: make user choose path
    lambda.path <- get_lambda_path(x = X, time = time, 
                                   status = status, 
                                   failcode = failcode, 
                                   alpha = alpha,
                                   m = 100)
    
    # fit competing risk model with elastic net penalty along the lambda path
    model.obj <- fastcmprsk::fastCrrp(fastcmprsk::Crisk(time, status, failcode = failcode) ~.,
                                      data = data.frame(time, status, X),
                                      penalty = "ENET",
                                      alpha = alpha,
                                      lambda = lambda.path,...)
    
    # get unexported function 'AIC.fcrrp'
    fnc <- utils::getFromNamespace("AIC.fcrrp", "fastcmprsk")
    
    # choose the penalty that minimizes AIC
    s <- unname(which.min(fnc(model.obj)))
    
    # store it
    lambda.min <- lambda.path[s]
    
    # store coefficient matrix as sparse matrix (for consistency with glmnet)
    coefs <- Matrix::Matrix(model.obj$coef[,s,drop=FALSE], sparse = TRUE, 
                            dimnames = list(colnames(X), "Estimate"))
    
    # get indices of kept variables (account for zero indexing)
    kept.vars     <- coefs@i + 1
    
    # get linear predictor 
    lp  <- as.numeric(X[,kept.vars,drop = FALSE] %*% coefs[kept.vars])
    
  } else{
    
    # fit competing risk model without penalty
    model.obj <- cmprsk::crr(ftime = time, fstatus = status, 
                             cov1 = X, failcode = failcode,...)
    
    # store coefficient matrix as sparse matrix (for consistency with glmnet)
    coefs <- Matrix::Matrix(model.obj$coef, sparse = TRUE, 
                            dimnames = list(colnames(X), "Estimate"))
    
    # get linear predictor 
    lp  <- as.numeric(X %*% coefs)
    
    # empty lambda (no regularization)
    lambda.min <- NULL
    
  } # IF
    
  
  # get the survival functions    
  surv.obj <- survival_cmprsk(time = time, status = status, lp = lp, 
                  prep_predict_object = NULL,
                  failcode = failcode)
  
  # return
  return(structure(list(
    risk = 1.0 - surv.obj$surv(time_eval = time_eval),
    linear_predictor = lp,
    coefficients = coefs,
    model = model.obj,
    lambda_min = lambda.min,
    funs = surv.obj
  ), class = "baseline_survival"))
  
} # FUN


baseline_survival_nocmprsk <- function(X, 
                                       status,
                                       time,
                                       time_eval = max(time),
                                       alpha = 1,
                                       failcode = 1, # dead argument here
                                       ...)
{
  
  if(!is.null(alpha)){
    
    # prepare dependent variable
    y <- survival::Surv(time = time, event = status)
    
    # fit models
    model.obj <- glmnet::cv.glmnet(x = X, y = y,
                                   family = "cox",
                                   type.measure = "deviance",
                                   alpha = alpha,...) 
    
    # get the coefficients
    coefs <- glmnet::coef.glmnet(model.obj, s = "lambda.min")
    colnames(coefs) <- "Estimate"
    
    # get indices of kept variables (account for zero indexing)
    kept.vars     <- coefs@i + 1
    
    # get linear predictor 
    lp  <- as.numeric(X[,kept.vars,drop = FALSE] %*% coefs[kept.vars])
    
    # get minimizing lambda
    lambda.min <- model.obj$lambda.min
  
  } else{
    
    # fit Cox PH model
    model.obj <- survival::coxph(survival::Surv(time = time, event = status)~.,
                                 data = data.frame(time, status, X), ...)
    
    # store coefficient matrix as sparse matrix (for consistency with glmnet)
    coefs <- Matrix::Matrix(model.obj$coefficients, sparse = TRUE, 
                            dimnames = list(colnames(X), "Estimate"))
    
    # get linear predictor 
    lp  <- as.numeric(X %*% coefs)
    
    # empty lambda parameter (no regularization)
    lambda.min <- NULL
    
  } # IF
  
  # get the survival functions    
  surv.obj <- survival(time = time, status = status, lp = lp, center = FALSE)
  
  # return
  return(structure(list(
    risk = 1.0 - surv.obj$surv(time_eval = time_eval),
    linear_predictor = lp,
    coefficients = coefs,
    model = model.obj,
    lambda_min = lambda.min,
    funs = surv.obj
  ), class = "baseline_survival"))
  
} # FUN


#' Predict method for a \code{baseline_crss} object
#' 
#' @param object A \code{baseline_crss} object.
#' @param newX A numeric matrix at which predictions should be performed
#' @param shrunk Shall 1st (TRUE) or 2nd (FALSE) stage be used for predictions?
#' @param ... Additional parameters to be passed down
#' 
#' @return A matrix of risk predictions
#' 
#' @export
predict.baseline_risk <- function(object, newX, shrunk = FALSE, ...)
{
  ## input checks
  if(!inherits(x = object, what = "baseline_crss", which = FALSE))
  {
    stop("object must be an instance of baseline_crss()")
  }
  
  InputChecks_newX(newX)
  newX <- check_and_adjust_newX(newX = newX,
                                covariates = colnames(object$inputs$X))
  
  ## predict
  predict_baseline_crss_NoChecks(object = object, newX = newX, shrunk = shrunk, ... = ...)
  
} # FUN


predict_baseline_crss_NoChecks <- function(object, newX, shrunk, ...)
{
  
  mod_type <- ifelse(shrunk, "full", "reduced")
  retained <- rownames(object$coefficients$reduced) # gives us coefs with nonzero 1st stage estimates
  
  if("(Intercept)" %in% retained)
  {
    retained0 <- retained[-1L]
    intercept <- object$coefficients[[mod_type]]["(Intercept)", "Estimate"]
  } else{
    intercept <- 0.0
  }
  
  if(identical(length(retained0), 0L))
  {
    # case 1: no variables were retained
    out <- rep(intercept, nrow(newX))
  } else
  {
    # case 2: at least one variable was retained
    cf <- as.numeric(object$coefficients[[mod_type]][retained, "Estimate"])
    lp <- cbind(1.0, newX[,retained0,drop = FALSE]) %*% cf
    out <- as.numeric(stats::plogis(lp))
    
  } # IF
  
  return(matrix(out))

} # FUN


predict_baseline_survival_NoChecks <- function(object, newX, time_eval, ...)
{
  ### get the linear predictor (lp) at 'newX' 
  ## case 1: glmnet object
  if(inherits(x = object$model, what = "cv.glmnet"))
  {
    lp <- glmnet:::predict.cv.glmnet(
      object = object$model,
      newx = newX, s = "lambda.min", 
      type = "link", ... = ...) 
    
  } else if(inherits(x = object$model, what = c("coxph", "crr"))){
    
    ## case 2: no regularization
    # applies to 'coxph' and 'crr' objects
    lp <- newX %*% as.numeric(object$coefficients)
    
  } else{
    
    ## case 3: regularized competing risks model (class 'fcrrp')
    # get indices of kept variables (account for zero indexing)
    kept_vars <- object$coefficients@i + 1
    lp        <- newX[,kept_vars,drop = FALSE] %*% object$coefficients[kept_vars]
    
  } # IF
  
  ## convert lp to vector
  lp <- as.numeric(lp)
  
  ## calculate baseline survival probability at time_eval
  s0 <- object$funs$basesurv(time_eval = time_eval)
  
  ## calculate risk as 1 - {survival probability}
  risk <- 1.0 - s0 * exp(lp)
  
  return(matrix(risk))
  
} # FUN


#' Predict method for a \code{baseline_survival} object
#' 
#' @param object A \code{baseline_survival} object.
#' @param newX A numeric matrix at which predictions should be performed.
#' @param time_eval The time at which baseline risk shall be predicted. Must be nonnegative numeric vector of length one.
#' @param ... Additional parameters to be passed down
#' 
#' @return A matrix of risk predictions
#' 
#' @export
predict.baseline_survival <- function(object, newX, time_eval, ...)
{
  ## input checks
  if(!inherits(x = object, what = "baseline_survival", which = FALSE))
  {
    stop("object must be an instance of baseline_survival()")
  }
  
  InputChecks_newX(newX)
  InputChecks_newX_X(newX = newX, object = object, survival = TRUE)
  stopifnot(length(time_eval) == 1L & is.numeric(time_eval))
  
  ## predict
  predict_baseline_survival_NoChecks(object = object, 
                                     newX = newX,
                                     time_eval = time_eval,
                                     ... = ...)
  
} # FUN