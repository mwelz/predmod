#' Performs risk modeling by logistic regression.
#' 
#' The logistic regression model is specified as
#' \deqn{logit(Pr[y=1|z,w]) = \alpha + \beta_1 w + \beta_2 z + \beta_3 z * w ,}
#' where \eqn{z} is the typically equal to the linear predictor of a baseline risk model. If \code{constant = TRUE}, it is enforced that \eqn{\beta_3 = 0}. 
#' 
#' @param X Matrix of fixed covariates.
#' @param status Numeric vector with a unique code for each failure type. Code 0 denotes a censored observation.
#' @param w Binary vector of treatment assignment status. Equal to 1 for treatment group and 0 for control group.
#' @param z The \code{z} to be used in the logistic regression model above. If \code{NULL} (default), then the linear predictor of a baseline risk model is used as \code{z}.
#' @param alpha The elasticnet mixing parameter for regularization in a potential baseline risk model to obtain the linear predictor to be used as \code{z}. Only applicable if \code{z = NULL}. See \code{\link{baseline_risk}} for details.
#' @param failcode Code of status that denotes the failure type of interest. Default is one.
#' @param constant Shall \eqn{\beta_3 = 0} be enforced in the logistic regression model above? If \code{TRUE}, then this is equivalent to assuming that there is a constant treatment effect. Default is \code{FALSE}.
#' @param LRT shall likelihood ratio test performed to test if there is a constant treatment effect? Default is \code{FALSE}.
#' @param significance_level Significance level of possible likelihood ratio test. Only applicable if \code{LRT = TRUE}.
#' @param glm_data TODO
#' @param ... Additional arguments to be passed.
#' 
#' @return A \code{predmod_ordinary} object.
#' 
#' @export
risk_model <- function(X, 
                       status, 
                       w, 
                       z = NULL, 
                       alpha = 1,
                       failcode = 1,
                       constant = FALSE,
                       LRT = TRUE,
                       significance_level = 0.05,
                       glm_data = FALSE,
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
  
  ## stage 1
  if(is.null(z)){
    
    # in case no z is supplied, fit a baseline risk model
    stage1 <- baseline_risk(X = X, status = status,
                            alpha = alpha, glm_data = glm_data,
                            failcode = failcode,...)
    
    # take the baseline model's linear predictor as z
    z <- stage1$linear_predictor
    
    # extract baseline risk
    baseline.risk <- stage1$risk
    
  } else{
    
    # in case values for z are supplied, don't fit baseline risk model
    stage1 <- baseline.risk <- NULL
    
  } # IF
  
  
  ## stage 2
  if(LRT)
  {
    fit <- LRT_crss(status = status_bin, 
                    w = w, 
                    w_flipped = ifelse(w == 1, 0, 1), 
                    z = z, 
                    significance_level = significance_level,
                    glm_data = glm_data,
                    ... = ...)
    
    # extract risk estimates
    risk_reg <- fit$models$accepted$risk$regular
    risk_rev <- fit$models$accepted$risk$counterfactual
    
    # organize output
    stage2 <- list(
      models = list(accepted = fit$models$accepted$model, 
                    rejected = fit$models$rejected$model),
      deviance = fit$deviance,
      pval_LRT = fit$pval_LRT,
      decision = fit$decision
    )
    
  } else{
    fit <- risk_model_stage2(status_bin = status_bin, 
                                w = w, 
                                w_flipped = ifelse(w == 1, 0, 1), 
                                z = z, 
                                constant = constant,
                                glm_data = glm_data,
                                ...  = ...)
    
    # extract risk estimates
    risk_reg <- fit$risk$regular
    risk_rev <- fit$risk$counterfactual
    
    # organize output
    stage2 <- list(
      models = list(accepted = fit$model, rejected = NULL),
      deviance = NULL,
      pval_LRT = NULL,
      decision = NULL
    )
    
  } # IF
  
  
  # calculate predicted benefits
  benefits <- get_predicted_benefits(risk_reg = risk_reg, 
                                     risk_rev = risk_rev,
                                     w = w)
  
  # return
  return(structure(list(benefits = benefits,
                        coefficients = list(baseline = stage1$coefficients,
                                            stage2 = list(
                                              accepted = get_coefs(stage2$models$accepted),
                                              rejected = get_coefs(stage2$models$rejected)
                                            )),
                        risk = list(baseline = baseline.risk,
                                    regular = as.matrix(risk_reg),
                                    counterfactual = as.matrix(risk_rev)),
                        models = list(baseline = stage1, 
                                      stage2 = stage2),
                        inputs = list(status = status, status_bin = status_bin,
                                      w = w, failcode = failcode, z = z, 
                                      constant = constant, alpha = alpha, 
                                      covariates = colnames(X))
  ), 
  class = "risk_model_crss"))
  
} # FUN


#' (For internal use only.) Fits the second stage risk model. 
#' 
#' @param status_bin A \strong{binary} vector of status. Zero is survivor, one is failure.
#' @param A binary treatment assignment status.
#' @param w_flipped Treatment assignment status, but flipped .
#' @param constant Do we assume constant treatment effect?
#' @param glm_data TODO
#' @param ... Additional arguments,
#' 
#' @return A list with a \code{glm} object as well as the two risk predictions.
#'
#' @noRd
risk_model_stage2 <- function(status_bin, w, z,
                              w_flipped = ifelse(w == 1, 0, 1), 
                              constant = FALSE,
                              glm_data = FALSE,
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
  model <- stats::glm(status_bin~., 
                      family = stats::binomial(link = "logit"), 
                      data = as.data.frame(X_stage2, status_bin),
                      x = FALSE, y = FALSE, model = FALSE, ...)
  
  # get the responses with the regular w ( = F_logistic(x'beta + z))
  risk_reg <- as.numeric(stats::plogis(cbind(1.0, X_stage2) %*% model$coefficients))
  
  # get responses with flipped W
  risk_rev <- as.numeric(stats::plogis(cbind(1.0, X_stage2_rev) %*% model$coefficients))
  
  # drop data from glm object if desired
  if(!glm_data) model$data <- NULL
  
  # return
  return(list(model = model,
              risk = list(regular = risk_reg, counterfactual = risk_rev)
  ))
} # FUN
