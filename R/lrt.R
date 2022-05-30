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
