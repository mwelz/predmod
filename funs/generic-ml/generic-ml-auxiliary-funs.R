#' Estimates the BLP parameters based on the main sample M. 
#' 
#' @param Z an ( _|M|_ x _d_) matrix or data frame of covariates
#' @param D a binary vector of treatment status of length _|M|_
#' @param Y a vector of responses of length _|M|_
#' @param propensity.scores a vector of propensity scores of length _|M|_
#' @return BLP coefficients with inference statements
#' 
#' @export
#' 
#' TODO: implement same with HT transformation! 
get.BLP.params.classic <- function(Z, D, Y, propensity.scores){
  
  # prepare weights
  weights <- 1 / (propensity.scores * (1 - propensity.scores))
  
  # prepare covariate matrix
  X <- data.frame(B = proxy.baseline, 
                  S = proxy.cate,
                  beta1 = D - propensity.scores, 
                  beta2 = (D - propensity.scores) * (proxy.cate - mean(proxy.cate))) 
  
  # fit weighted linear regression by OLS
  blp.obj <- lm(Y ~., data = data.frame(Y, X), weights = weights)
  
  # extract coefficients
  coefficients     <- summary(blp.obj)$coefficients
  
  # inference on beta2: test the null that it is 1) = 0, 2) = 1.
  beta2.inference <- matrix(NA_real_, 2, 2)
  rownames(beta2.inference) <- c("H0: beta2 = 0", "H0: beta2 = 1")
  colnames(beta2.inference) <- c("t value", "Pr(>|t|)")
  beta2.inference[1,] <- coefficients["beta2", c("t value", "Pr(>|t|)")]
  beta2.inference[2, "t value"] <- 
    (coefficients["beta2", "Estimate"] - 1) / coefficients["beta2", "Std. Error"]  
  beta2.inference[2, "Pr(>|t|)"] <- 
    2 * pt(beta2.inference[2, "t value"], df = blp.obj$df.residual, lower.tail = FALSE)
  
  return(list(lm.obj = blp.obj, 
              blp.coefficients = blp.obj$coefficients[c("beta1", "beta2")],
              coefficients = coefficients,
              beta2.inference = beta2.inference))
  
} # END FUN