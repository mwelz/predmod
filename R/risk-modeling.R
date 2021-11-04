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
                            alpha = alpha, 
                            failcode = failcode,...)
    
    # take the baseline model's linear predictor as z
    z <- stage1$linear_predictor
    
    # estimate concordance on 1st stage
    C_outcome_stage1 <- C_outcome(y = status_bin, stage1$risk)
    
    # extract baseline risk
    baseline.risk <- stage1$risk
    
  } else{
    
    # in case values for z are supplied, don't fit baseline risk model
    stage1 <- C_outcome_stage1 <- baseline.risk <- NULL
    
  } # IF
  

  ## stage 2
  stage2 <- risk_model_stage2(status_bin = status_bin, 
                              w = w, 
                              w_flipped = ifelse(w == 1, 0, 1), 
                              z = z, 
                              constant = constant)

  # extract risk estimates
  risk_reg <- stage2$risk$regular
  risk_rev <- stage2$risk$counterfactual
  
  # calculate predicted benefits
  benefits <- get_predicted_benefits(risk_reg = risk_reg, 
                                     risk_rev = risk_rev,
                                     w = w)

  # return
  return(structure(list(benefits = benefits,
                        coefficients = summary(stage2$model)$coefficients,
                        risk = list(baseline = baseline.risk,
                                    regular = risk_reg,
                                    counterfactual = risk_rev),
                        concordance = list(outcome = C_outcome_stage1,
                                           benefit = C_benefit(y = status_bin, 
                                                               w = w,
                                                               pred_ben = benefits$absolute)),
                        models = list(stage1 = stage1, stage2 = stage2),
                        inputs = list(status = status, status_bin = status_bin,
                                      w = w, failcode = failcode, z = z, 
                                      constant = constant, alpha = alpha)
                        ), 
                   class = "predmod_ordinary"))
  
} # FUN





#' (For internal use only.) Fits the second stage risk model. 
#' 
#' @param status_bin A \strong{binary} vector of status. Zero is survivor, one is failure.
#' @param A binary treatment assignment status.
#' @param w_flipped Treatment assignment status, but flipped .
#' @param constant Do we assume constant treatment effect?
#' @param ... Additional arguments,
#' 
#' @return A list with a \code{glm} object as well as the two risk predictions.
#'
#' @noRd
risk_model_stage2 <- function(status_bin, w, z,
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
  model <- stats::glm(status_bin~., 
                      family = stats::binomial(link = "logit"), 
                      data = as.data.frame(X_stage2, status_bin),...)
  
  # get the responses with the regular w ( = F_logistic(x'beta + z))
  risk_reg <- as.numeric(stats::predict.glm(model, type = "response"))
  
  # get responses with flipped W
  risk_rev <- as.numeric(stats::plogis(cbind(1, X_stage2_rev) %*% model$coefficients))
  
  # return
  return(list(model = model,
              risk = list(regular = risk_reg, counterfactual = risk_rev)
              ))
} # FUN