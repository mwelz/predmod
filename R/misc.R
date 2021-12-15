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
#' @param x A predmod object
#' @param subset The indices of the subgroup of interest
#' @param relative Shall relative ATE be calculated?
#' @param benefits_risk Logical. If \code{TRUE}, then the failure-risk-based benefits are used (only applicable to survival models). Default is \code{FALSE}.
#' @param time_eval Only applicable if \code{benefits_risk = TRUE}. Time at which to evaluate the failure risk predictions.
#' 
#' @export
average_treatment_effect <- function(x, 
                                     subset = NULL, 
                                     relative = FALSE,
                                     benefits_risk = FALSE,
                                     time_eval = NULL){
  
  clss <- class(x)
  w    <- x$inputs$w
  stopifnot(clss %in% c("predmod_ordinary", "predmod_survival"))
  
  # prepare subset object
  if(!is.null(subset)){
    stopifnot(is.numeric(subset))
  } else {
    subset <- 1:length(w)
  } # IF
  
  # prepare time_eval object
  if(is.null(time_eval)){
    
    # if no value provided, take time_eval used in model fitting
    if(clss == "predmod_survival"){
      time_eval <- x$input$time_eval
    } # IF clss
    
  } else {
    
    stopifnot(length(time_eval) == 1L)
    
  } # IF time_eval
  
  
  if(!relative){
    
    # in case of absolute effect, we need to decompose the effect
    # so that we can apply a two-sample test
    
    if(clss == "predmod_ordinary"){
      
      # case 1: predmod_ordinary
      x_reg <- x$risk$regular
      x_rev <- x$risk$counterfactual
      
    } else{
      
      # case 2: predmod_survival
      if(benefits_risk){
        
        # case 2.1: benefits concern risk at time of interest
        x_reg <- 1.0 - x$funs$regular$surv(time_eval)
        x_rev <- 1.0 - x$funs$counterfactual$surv(time_eval)
        
      } else{
        
        # case 2.2: benefits concern failure times
        x_reg <- x$failure$regular
        x_rev <- x$failure$counterfactual
        
      } # IF benefits_risk
    } # IF class
    
    
    # adjust signs to ensure that x_reg - x_rev = (predicted absolute benefit) 
    x_reg[w == 0] <- -x_reg[w == 0]
    x_rev[w == 0] <- -x_rev[w == 0]
    
    # perform t-test
    t <- stats::t.test(x = x_reg[subset], y = x_rev[subset], paired = FALSE, var.equal = FALSE)
    
    # get ATE
    ate <- unname(t$estimate[1] - t$estimate[2])
    sderr <- t$stderr
    
  } else {
    
    # relative effect: here we can simply use the direct estimated benefits
    if(clss == "predmod_ordinary"){
      
      # case 1: benefits concern risk in cross-sectional model
      ate <- mean(x$benefits$relative[subset])
      
    } else{
      
      if(benefits_risk){
        
        # case 2.1: benefits concern risk at time of interest
        x_reg      <- 1.0 - x$funs$regular$surv(time_eval)
        x_rev      <- 1.0 - x$funs$counterfactual$surv(time_eval)
        rr         <- x_reg / x_rev
        rr[w == 0] <- 1 / rr[w == 0]
        ate        <- mean(rr[subset])
        
      } else{
        
        # case 2.2: benefits concern failure time
        ate <- mean(x$benefits$relative[subset])
        
      } # IF benefits_risk
    } # IF clss
    
    # no SE can be computed for relative risk TODO: maybe via bootstrap
    sderr <- NA_real_
    
  } # IF !relative
 
  # return
  return(structure(c(ate, sderr), names = c("ATE", "Std. Error")))

} # FUN


#' Partition a vector into quantile groups
#'
#' Partitions a vector into quantile groups and returns a logical matrix indicating group membership.
#'
#' @param x A numeric vector to be partitioned.
#' @param cutoffs A numeric vector denoting the quantile cutoffs for the partition. Default are the quartiles: \code{c(0.25, 0.5, 0.75)}.
#'
#' @return
#' An object of type \code{quantile_group}, which is a logical matrix indicating group membership.
#'
#' @examples
#' set.seed(1)
#' x <- runif(100)
#' cutoffs <- c(0.25, 0.5, 0.75)
#' quantile_group(x, cutoffs)
#'
#' @export
quantile_group <- function(x,
                           cutoffs = c(0.25, 0.5, 0.75)){
  
  
  # input checks
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(cutoffs))
  stopifnot(0 < min(cutoffs) & max(cutoffs) < 1)
  
  # return
  quantile_group_NoChecks(x = x, cutoffs = cutoffs)
  
  
} # FUN


# same as above, just w/o input checks
quantile_group_NoChecks <- function(x = x,
                                    cutoffs = cutoffs){
  
  # get quantiles
  q         <- stats::quantile(x, cutoffs)
  q         <- c(-Inf, q, Inf)
  
  # check if breaks are unique: if x exhibits low variation, there might be empty quantile bins, which can cause an error in the cut() function. In this case, we add random noise to x to induce variation. NB: this bug has been spotted and fixed by Lucas Kitzmueller. All credits for this fix go to him!
  if(length(unique(q)) != length(q)){
    # specify standard deviation of the noise (x may have zero variation)
    sd <- ifelse(stats::var(x) == 0, 0.001, sqrt(stats::var(x) / 20))
    # add noise and updare quantiles
    x <- x + stats::rnorm(length(x), mean = 0, sd = sd)
    q <- stats::quantile(x, cutoffs)
    q <- c(-Inf, q, Inf)
  } # IF
  
  groups    <- as.character(cut(x, breaks = q, include.lowest = TRUE, right = FALSE, dig.lab = 3))
  group.nam <- unique(groups)
  group.nam <- group.nam[order(
    as.numeric(substr(sub("\\,.*", "", group.nam), 2, stop = 1e8L)),
    decreasing = FALSE)] # ensure the order is correct
  
  # get the grouping matrix
  group.mat <- sapply(1:length(group.nam), function(j) groups == group.nam[j])
  colnames(group.mat) <- gsub(",", ", ", gsub(" ", "", group.nam))
  
  # return
  return(structure(group.mat, type = "quantile_group"))
  
} # FUN

intervals.quantile <- function(cutoffs){
  K   <- length(cutoffs) + 1
  nam <- rep(NA_character_, K)
  
  for(j in 1:K){
    if(j == 1){
      nam[j] <- paste0("<", 100*cutoffs[j], "%")
    } else if (j == K){
      nam[j] <- paste0(">=", 100*cutoffs[j-1], "%")
    } else{
      nam[j] <- paste0("[", 100*cutoffs[j-1], ",", 100*cutoffs[j], ")%")
    }
  } # FOR
  
  nam
  
} # FOR
