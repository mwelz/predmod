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


# Taken from https://github.com/grf-labs/grf/blob/bf691429b5714385da5f41f6b5db1525184b300d/experiments/csf/comparison_estimators.R#L8 , credits to to the authors of the grf package!
expected_survival <- function(S.hat, Y.grid) {
  grid.diff <- diff(c(0, Y.grid, max(Y.grid)))
  
  c(cbind(1, S.hat) %*% grid.diff)
}


#' computes the average treatment effect (ATE) of a predmod object
#' 
#' @param predmod A predmod object
#' @param subset The indices of the subgroup of interest
#' 
#' @export
average_treatment_effect <- function(predmod, subset = NULL){
  
  if(!is.null(subset)) stopifnot(is.numeric(subset))
  clss <- class(predmod)
  
  if(clss == "predmod_survival"){
    
    x_reg <- predmod$failure$regular
    x_rev <- predmod$failure$counterfactual
    
  } else if(clss == "predmod_ordinary"){
    
    x_reg <- predmod$risk$regular
    x_rev <- predmod$risk$counterfactual
    
  } else stop("This function is only impplemented for classes 'predmod_ordinary' and 'predmod_survival'.")
  
  w <- predmod$inputs$w
  
  if(is.null(subset)){
    idx <- 1:length(w)
  } else{
    idx <- subset
  }
  
  # adjust signs to ensure that x_reg - x_rev = (predicted absolute benefit) 
  x_reg[w == 0] <- -x_reg[w == 0]
  x_rev[w == 0] <- -x_rev[w == 0]
  
  # perform t-test
  t <- stats::t.test(x = x_reg[idx], y = x_rev[idx], paired = FALSE, var.equal = FALSE)
  ate <- unname(t$estimate[1] - t$estimate[2])
  
  # return
  return(structure(c(ate, t$stderr), names = c("ATE", "Std. Error")))

} # FOR