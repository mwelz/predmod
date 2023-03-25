#' calculate likelihood ratio test statistic (LRT) for cross-sectional models
#' 
#' Performs a likelihood ratio test between models with and without constant treatment effect using the lmtest package
#' @param status A \strong{binary} vector of status. Zero is survivor, one is failure.
#' @param w A binary treatment assignment status.
#' @param w_flipped Treatment assignment status, but flipped.
#' @param z The linear predictor of the first stage model
#' @param significance_level Significance level of test
#' @param ... Additional arguments
#' @references  Zeileis A, Hothorn T (2002). “Diagnostic Checking in Regression Relationships.” R News, 2(3), 7–10. https://CRAN.R-project.org/doc/Rnews/. 
#' @noRd
LRT_crss <- function(status, w, w_flipped, z, significance_level, glm_data, ...)
{
  
  ## fit restricted and full model
  stage2_0 <- 
    risk_model_stage2(status_bin = status, 
                      w = w, 
                      w_flipped = w_flipped, 
                      z = z, 
                      constant = TRUE, 
                      glm_data = glm_data,
                      ... = ...)
  
  stage2_full <- 
    risk_model_stage2(status_bin = status, 
                      w = w, 
                      w_flipped = w_flipped, 
                      z = z, 
                      constant = FALSE, 
                      glm_data = glm_data,
                      ... = ...)
  
  mod_0    <- stage2_0$model
  mod_full <- stage2_full$model
  
  ## calculate deviance
  lmtest_obj <- lmtest::lrtest(mod_full, mod_0)
  deviance   <- 2.0 * (lmtest_obj$LogLik[1] - lmtest_obj$LogLik[2])
  
  ## return
  LRT_return(stage2_full = stage2_full,
             stage2_0 = stage2_0,
             deviance = deviance,
             significance_level = significance_level)
  
} # FUN


<<<<<<< HEAD
=======

>>>>>>> master
# helper function that formats the returns of the LRT function
LRT_return <- function(stage2_full, stage2_0, deviance, significance_level)
{
  
  # calculate p-value of the LRT test
  pval <- stats::pchisq(deviance, df = 1L, 
                        lower.tail = FALSE)
  
  if(pval < significance_level)
  {
    ## case 1: reject the null of the models having equal fit
    # in this case, go for the full model
    acc <- stage2_full
    rej <- stage2_0
    dec <- "full"
    
  } else{
    
    ## case 2: do not reject null of equal fit
    # in this case, go for smaller model
    acc <- stage2_0
    rej <- stage2_full
    dec <- "reduced"
  }
  
  # return both models, deviance, and p-value
  out <- list(
    models = list(accepted = acc, rejected = rej),
    deviance = deviance,
    pval_LRT = pval,
    decision = dec
  )
  
  return(out)
  
} # FUN