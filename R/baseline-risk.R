#' baseline risk prediction via penalized logistic regression (treatment assignment is not used)
#' 
#' Hello there
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
