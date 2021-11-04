# Benefits lines line 424: add
# 
# # add noise in case of low variation
# 
# if(var(rel.pred.ben) < 0.005){
#   
#   rel.pred.ben <- rel.pred.ben + rnorm(length(rel.pred.ben), 0, sqrt(var(y)/20))
#   
# }
# 
# if(var(abs.pred.ben) < 0.005){
#   
#   abs.pred.ben <- abs.pred.ben + rnorm(length(abs.pred.ben), 0, sqrt(var(y)/20))
#   
# }


#' (For internal use only.) Returns the lambda path for regularized Cox models as used in \href{https://cran.r-project.org/web/packages/glmnet/}{glmnet} package.
#' 
#' @param x A data matrix.
#' @param time Right-censored time at risk.
#' @param status Numeric vector with a unique code for each failure type and a separate code for censored observations.
#' @param failcode Code of \code{status} that denotes the failure type of interest.
#' @param alpha The elasticnet mixing parameter, with \eqn{0\le\alpha\le 1}.
#' The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty.
#' @param m Length of the path.
#' 
#' @noRd
get_lambda_path <- function(x, time, status, failcode = 1, alpha = 1, m = 100){
  
  # make status binary: one if failure due to cause of interest
  status.failcode <- ifelse(status == failcode, 1, 0)
  
  # get unexported function 'get_cox_lambda_max'
  fnc <- utils::getFromNamespace("get_cox_lambda_max", "glmnet")
  
  # get lambda.max by using standardized design matrix
  lambda.max <- fnc(x = scale(x, TRUE, TRUE),
                    y = survival::Surv(time, status.failcode),
                    alpha = alpha)
  
  # get the lambda path as suggested in Simon et al. (2011, JSS)
  eps <- ifelse(nrow(x) < ncol(x), 0.05, 0.0001)
  lambda.max * sapply(0:m, function(j) eps^(j/m))
  
} # FUN