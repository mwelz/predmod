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
  risk_reg <- 
    unname(
      stats::predict.glm(
        object = fits$models$reduced,
        newdata = as.data.frame(X_full[,fits$kept_vars, drop = FALSE]), type = "response"))
  
  risk_rev <- 
    unname(
      stats::predict.glm(
        object = fits$models$reduced,
        newdata = as.data.frame(X_full_rev[,fits$kept_vars, drop = FALSE]), type = "response"))

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
    br <- as.matrix(baseline_risk)
    mod_baseline <- NULL
  } # IF
  
  # return
  return(structure(list(benefits = benefits,
                        coefficients = list(baseline = mod_baseline$coefficients,
                                            full = fits$coefficients$full,
                                            reduced = summary(fits$models$reduced)$coefficients),
                        risk = list(baseline = br,
                                    regular = as.matrix(risk_reg),
                                    counterfactual = as.matrix(risk_rev)),
                        models = list(baseline = mod_baseline, 
                                      full = fits$models$full,
                                      reduced = fits$models$reduced),
                        inputs = list(status = status, status_bin = status_bin,
                                      w = w, failcode = failcode, alpha = alpha, 
                                      alpha_baseline = alpha_baseline, 
                                      interacted = interacted,
                                      covariates = colnames(X))
  ), 
  class = "effect_model_crss"))
  
  
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
                              family =  stats::binomial(link = "logit"), x = FALSE, y = FALSE)
  
  # store coefficient matrix as sparse matrix (for consistency with glmnet)
  coefs_reduced <- Matrix::Matrix(model_reduced$coefficients, sparse = TRUE, 
                                  dimnames = list(c("(Intercept)", colnames(X_reduced)), "Estimate"))
  
  return(list(models = list(full = model_full, reduced = model_reduced),
              coefficients = list(full = coefs_full, reduced = coefs_reduced),
              kept_vars = kept.vars,
              lambda_min = model_full$lambda.min))
  
} # FUN
