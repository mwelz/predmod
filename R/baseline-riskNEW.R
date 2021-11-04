#' baseline risk prediction via penalized logistic regression (treatment assignment is not used)
#' 
#' @param X Matrix of covariates
#' @param status Vector of binary responses
#' @param alpha The elasticnet mixing parameter for regularization, with \eqn{0\le\alpha\le 1}.
#' The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1} (default) is the lasso penalty, and \code{alpha=0} the ridge penalty. If \code{NULL}, no regularization is used.
#' @param ... Additional arguments to be passed
#' 
#' @export
baseline_risk <- function(X,
                          status, 
                          alpha = 1,
                          ...)
{
  
  # input checks
  InputChecks_NA(list(X, status))
  CheckInputs_X(X)
  InputChecks_equal.length2(X, status)
  InputChecks_Y_binary(status)
  
  # assign variable names if there are none
  if(is.null(colnames(X))) colnames(X) <- paste0("V", 1:ncol(X)) 
  
  if(is.null(alpha)){
    
    ## case 1: no regularization
    
    # fit the model
    model.obj <- stats::glm(status~., 
                            family = stats::binomial(link = "logit"), 
                            data = data.frame(status, X), ...)
    
    # store coefficient matrix as sparse matrix (for consistency with glmnet)
    coefs <- Matrix::Matrix(model.obj$coefficients, sparse = TRUE, 
                            dimnames = list(colnames(X), NULL))
    
    # get linear predictor 
    lp  <- as.numeric(cbind(1,X) %*% coefs)
    
    # empty lambda parameter (no regularization)
    lambda.min <- NULL
    
  } else{
    
    ## case 2: regularization
  
    # fit model
    model.obj     <- glmnet::cv.glmnet(X, status, 
                                       family = "binomial", 
                                       alpha = alpha, ...)
    
    # get coefficients at best lambda
    coefs <- glmnet::coef.glmnet(model.obj, s = "lambda.min")
    
    # get indices of kept variables (account for zero indexing)
    kept.vars     <- coefs@i + 1
    
    # get linear predictor 
    lp  <- as.numeric(cbind(1,X)[,kept.vars,drop = FALSE] %*% coefs[kept.vars])
    
    # get minimizing lambda
    lambda.min <- model.obj$lambda.min
    
  } # IF
  
  
  # return
  return(structure(list(
    risk = stats::plogis(lp),
    linear_predictor = lp,
    coefficients = coefs,
    model = model.obj,
    lambda = lambda.min
  ), class = "baseline_risk"))
  
} # FUN



#' baseline risk prediction via penalized logistic regression (treatment assignment is not used)
#' 
#' @param X Matrix of covariates
#' @param status Vector of binary responses
#' @param time Time at risk
#' @param time_eval timat at which survival shall be evaluated
#' @param alpha The elasticnet mixing parameter for regularization, with \eqn{0\le\alpha\le 1}.
#' The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1} (default) is the lasso penalty, and \code{alpha=0} the ridge penalty. If \code{NULL}, no regularization is used.
#' @param failcode Code of status that denotes the failure type of interest
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
  
  # TODO: fill in; this is a wrapper function
  # input checks
  InputChecks_NA(list(X, status))
  CheckInputs_X(X)
  InputChecks_equal.length2(X, status)
  InputChecks_Y_binary(status)
  stopifnot(length(time_eval) == 1)
  
  # assign variable names if there are none
  if(is.null(colnames(X))) colnames(X) <- paste0("V", 1:ncol(X)) 
  
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
                                      lambda = lambda.path)
    
    # choose the penalty that minimizes AIC
    s <- unname(which.min(fastcmprsk:::AIC.fcrrp(model.obj)))
    
    # store it
    lambda.min <- lambda.path[s]
    
    # store coefficient matrix as sparse matrix (for consistency with glmnet)
    coefs <- Matrix::Matrix(model.obj$coef[,s,drop=FALSE], sparse = TRUE, 
                            dimnames = list(colnames(X), NULL))
    
    # get indices of kept variables (account for zero indexing)
    kept.vars     <- coefs@i + 1
    
    # get linear predictor 
    lp  <- as.numeric(cbind(1,X)[,kept.vars,drop = FALSE] %*% coefs[kept.vars])
    
  } else{
    
    # fit competing risk model without penalty
    model.obj <- cmprsk::crr(ftime = time, fstatus = status, 
                             cov1 = X, failcode = failcode)
    
    # store coefficient matrix as sparse matrix (for consistency with glmnet)
    coefs <- Matrix::Matrix(model.obj$coef, sparse = TRUE, 
                            dimnames = list(colnames(X), NULL))
    
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
    lambda = lambda.min,
    funs = surv.obj
  ), class = "baseline_survival"))
  
} # FUN

