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
                       LRT = FALSE,
                       significance_level = 0.05,
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
    
    # extract baseline risk
    baseline.risk <- stage1$risk
    
  } else{
    
    # in case values for z are supplied, don't fit baseline risk model
    stage1 <- baseline.risk <- NULL
    
  } # IF
  
  
  ## stage 2
  if(LRT)
  {
    stage2 <- LRT(status = status_bin, 
                  w = w, 
                  w_flipped = ifelse(w == 1, 0, 1), 
                  z = z, 
                  significance_level = significance_level,
                  ... = ...)
  } else{
    stage2 <- risk_model_stage2(status_bin = status_bin, 
                                w = w, 
                                w_flipped = ifelse(w == 1, 0, 1), 
                                z = z, 
                                constant = constant,
                                ...  = ...)
  } # IF
  
  
  # extract risk estimates
  risk_reg <- stage2$risk$regular
  risk_rev <- stage2$risk$counterfactual
  
  # calculate predicted benefits
  benefits <- get_predicted_benefits(risk_reg = risk_reg, 
                                     risk_rev = risk_rev,
                                     w = w)
  
  # return
  return(structure(list(benefits = benefits,
                        coefficients = list(baseline = stage1$coefficients,
                                            stage2 = get_coefs(stage2$model)),
                        risk = list(baseline = baseline.risk,
                                    regular = as.matrix(risk_reg),
                                    counterfactual = as.matrix(risk_rev)),
                        models = list(baseline = stage1$model, stage2 = stage2$model),
                        inputs = list(status = status, status_bin = status_bin,
                                      w = w, failcode = failcode, z = z, 
                                      constant = constant, alpha = alpha)
  ), 
  class = "risk_model_crss"))
  
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
                      data = as.data.frame(X_stage2, status_bin),
                      x = FALSE, y = FALSE, ...)
  
  # get the responses with the regular w ( = F_logistic(x'beta + z))
  risk_reg <- as.numeric(stats::predict.glm(model, type = "response"))
  
  # get responses with flipped W
  risk_rev <- as.numeric(stats::plogis(cbind(1, X_stage2_rev) %*% model$coefficients))
  
  # return
  return(list(model = model,
              risk = list(regular = risk_reg, counterfactual = risk_rev),
              deviance = NULL, pval_LRT = NULL
  ))
} # FUN


#' calculate likelihood ratio test statistic (LRT)
#' 
#' LRT between models with and without constant treatment effect
#' @param status A \strong{binary} vector of status. Zero is survivor, one is failure.
#' @param A binary treatment assignment status.
#' @param w_flipped Treatment assignment status, but flipped .
#' @param z TODO
#' @param significance_level significance level of test
#' @param ... Additional arguments
#' 
#' @noRd
LRT <- function(status, w, w_flipped, z, significance_level = 0.05, ...)
{
  
  ## fit restricted and full model
  stage2_0 <- 
    risk_model_stage2(status_bin = status, 
                      w = w, 
                      w_flipped = w_flipped, 
                      z = z, 
                      constant = TRUE, 
                      ... = ...)
  
  stage2_full <- 
    risk_model_stage2(status_bin = status, 
                      w = w, 
                      w_flipped = w_flipped, 
                      z = z, 
                      constant = FALSE, 
                      ... = ...)
  
  mod_0    <- stage2_0$model
  mod_full <- stage2_full$model
  
  ## calculate deviance and p-value of the LRT test
  lmtest_obj <- lmtest::lrtest(mod_full, mod_0)
  deviance   <- 2.0 * abs(diff(lmtest_obj$LogLik))
  pval       <- stats::pchisq(deviance, df = 1L, 
                              lower.tail = FALSE)
  
  if(pval < significance_level)
  {
    ## case 1: reject the null of the models having equal fit
    # in this case, go for the full model
    out <- stage2_full
  } else{
    
    ## case 2: do not reject null of equal fit
    # in this case, go for smaller model
    out <- stage2_0
  }
  
  # add deviance and p-value
  out$deviance <- deviance
  out$pval_LRT <- pval
  
  return(out)
  
} # FUN
