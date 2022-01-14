#' (For internal use only.) Returns predicted benefits (absolute and relative) based on risk estimates.
#' 
#' @param risk_reg Risk estimates based on regular treatment assignment.
#' @param risk_rev Risk estimates based on reversed treatment assignment
#' @param w Binary treatment assignment.
#' 
#' @noRd
get_predicted_benefits <- function(risk_reg, risk_rev, w){
  
  # check for equal length
  InputChecks_equal.length3(risk_reg, risk_rev, w)
  
  # absolute predicted benefit
  absolute <- risk_reg - risk_rev
  
  # adjust for signs
  absolute[w == 0] <- -absolute[w == 0]
  
  # relative predicted benefit
  relative <- risk_reg / risk_rev
  
  # adjust for signs
  relative[w == 0] <- 1 / relative[w == 0]
  
  # return
  return(list(absolute = absolute, relative = relative))
  
} # FUN


#' calculates the absolute observed benefit, which corresponds to the difference in \code{mean(y[W=w])}, and the associated confidence interval.
#' 
#' @param status vector of binary outcomes
#' @param w vector of binary treatment status
#' @param significance_level the significance level. Default is 0.05.
#'
#' @export
observed_benefit_absolute <- function(status, w, significance_level = 0.05){
  
  if(stats::var(status) == 0){
    status <- status + stats::rnorm(length(status), 0, sqrt(0.001))
    warning("There is zero variation in 'status', so we added Gaussian noise with variance 0.001")
  } # IF
  
  ttest  <- stats::t.test(status[w==1], status[w==0], 
                          paired = FALSE, 
                          var.equal = FALSE) 
  aob    <- unname(ttest$estimate[1] - ttest$estimate[2])
  se.aob <- unname(ttest$stderr) 
  
  # quantile of the standard normal distribution
  z <- stats::qnorm(1 - significance_level/2)
  
  return(c(estimate = aob, 
           ci_lower = aob - z * se.aob,
           ci_upper = aob + z * se.aob,
           stderr   = se.aob))
} # FUN


#' calculates the relative observed benefit, which corresponds to a risk ratio, and the associated confidence interval
#' 
#' @param status vector of binary outcomes
#' @param w vector of binary treatment status
#' @param significance_level the significance level. Default is 0.05.
#'
#' @export
observed_benefit_relative <- function(status, w, significance_level = 0.05){
  
  # contingency table with the following structure
  #      Y=1 Y=0
  #  W=1  a   b
  #  W=0  c   d
  a <- sum(w == 1 & status == 1)
  b <- sum(w == 1 & status == 0)
  c <- sum(w == 0 & status == 1)
  d <- sum(w == 0 & status == 0)
  
  # estimates for the probabilities
  p1.hat <- a / (a + b) # estimates Pr(Y = 1 | W = 1)
  p0.hat <- c / (c + d) # estimates Pr(Y = 1 | W = 0)
  
  # empirical risk ratio (the relative observed benefit) and its standard error
  rob    <- p1.hat / p0.hat
  se.rob <- sqrt(b / (a * (a + b)) + d / (c * (c + d)))
  
  # quantile of the standard normal distribution
  z <- stats::qnorm(1 - significance_level/2)
  
  return(c(estimate = rob, 
           ci_lower = rob * exp(-z * se.rob),
           ci_upper = rob * exp(z * se.rob),
           stderr   = se.rob))
} # FUN



#' calculates the odds ratio, and the associated confidence interval. 
#' 
#' @param status vector of binary outcomes
#' @param w vector of binary treatment status
#' @param significance_level the significance level. Default is 0.05.
#'
#' @export
odds_ratio <- function(status, w, significance_level = 0.05){
  
  # contingency table with the following structure
  #      Y=1 Y=0
  #  W=1  a   b
  #  W=0  c   d
  a <- sum(w == 1 & status == 1)
  b <- sum(w == 1 & status == 0)
  c <- sum(w == 0 & status == 1)
  d <- sum(w == 0 & status == 0)
  
  if(any(a == 0 | b == 0 | c == 0 | d == 0)){
    stop("There are zero-valued entries in the contingency table")
  }
  
  # estimates for the probabilities
  p1.hat <- a / (a + b) # estimates Pr(Y = 1 | W = 1)
  p0.hat <- c / (c + d) # estimates Pr(Y = 1 | W = 0)
  
  # empirical odds ratio and its standard error
  or <- (p1.hat / (1 - p1.hat) ) / (p0.hat / (1 - p0.hat))
  se.or <- sqrt(1/a + 1/b + 1/c + 1/d)
  
  # quantile of the standard normal distribution
  z <- stats::qnorm(1 - significance_level/2)
  
  return(c(estimate = or, 
           ci_lower = or * exp(-z * se.or),
           ci_upper = or * exp(z * se.or),
           stderr   = se.or))
} # FUN



#' performs inference on the per-individual predicted benefits
#' 
#' @param x A predmod object
#' @param subset The indices of the subgroup of interest
#' @param relative Shall relative ATE be calculated?
#' @param benefits_risk Logical. If \code{TRUE}, then the failure-risk-based benefits are used (only applicable to survival models). Default is \code{FALSE}.
#' @param time_eval Only applicable if \code{benefits_risk = TRUE}. Time at which to evaluate the failure risk predictions.
#' @param significance_level The significance level, default is 0.05
#' 
#' @export
predicted_benefit_inference <- function(x, 
                                        subset = NULL, 
                                        relative = FALSE,
                                        benefits_risk = FALSE,
                                        time_eval = NULL,
                                        significance_level = 0.05){
  
  ate_obj <- average_treatment_effect(x = x, 
                                      subset = subset,
                                      relative = relative,
                                      benefits_risk = benefits_risk, 
                                      time_eval = time_eval)
  
  se  <- unname(ate_obj["Std. Error"])
  ate <- unname(ate_obj["ATE"])

  # quantile of the standard normal distribution
  z <- stats::qnorm(1 - significance_level/2)
  
  return(c(estimate = ate, 
           ci_lower = ate - z * se,
           ci_upper = ate + z * se,
           stderr   = se))
} # FUN


#' Calculates group-level benefits from a prediction model, and the associated confidence intervals.
#' The returned benefits are the observed and predicted relative and absolute benefits as well as the odds ratio
#' 
#' @param x prediction model object
#' @param cutoffs the quantile cutoff points. Default is c(0.25, 0.5, 0.75), which yields the quartiles.
#' @param baseline_risk The baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{x} is used.
#' @param benefits_risk Logical. If \code{TRUE}, then the risk-based benefits are used (only applicable to survival models). Default is \code{FALSE}.
#' @param time_eval Only applicable if \code{benefits_risk = TRUE}. Time at which we evaluate the risk predictions.
#' @param significance_level the significance level. Default is 0.05.
#' 
#' @export
get_benefits <- function(x, 
                         cutoffs = c(0.25, 0.5, 0.75),
                         baseline_risk = NULL,
                         benefits_risk = FALSE,
                         time_eval = NULL,
                         significance_level = 0.05){
  
  # extract outcome and treatment status
  status <- x$inputs$status_bin
  w      <- x$inputs$w
  
  # specify baseline risk that shall be used for grouping
  if(is.null(baseline_risk)){
    baseline_risk <- as.numeric(x$risk$baseline)
  } # IF
  
  # group observations by their quantile of predicted baseline risk
  quantile.groups <- quantile_group_NoChecks(baseline_risk, cutoffs)
  
  ## calculate observed benefit and predicted benefit for each quantile group
  # initialize
  rel.obs.ben.mat           <- matrix(NA_real_, ncol(quantile.groups), 4)
  colnames(rel.obs.ben.mat) <- c("estimate", "ci_lower", "ci_upper", "stderr")
  rownames(rel.obs.ben.mat) <- colnames(quantile.groups)
  abs.obs.ben.mat  <- rel.obs.ben.mat
  abs.pred.ben.mat <- rel.obs.ben.mat
  rel.pred.ben.mat <- rel.obs.ben.mat
  or.mat           <- rel.obs.ben.mat
  
  for(i in 1:ncol(quantile.groups)){
    
    group <- which(quantile.groups[,i])
    
    # absolute observed benefit
    abs.obs.ben.mat[i, ] <- 
      observed_benefit_absolute(status[group], w[group], significance_level = significance_level)
    
    # absolute predicted benefit
    abs.pred.ben.mat[i, ] <- 
      predicted_benefit_inference(x = x, subset = group, 
                                  relative = FALSE, 
                                  benefits_risk = benefits_risk, 
                                  time_eval = time_eval, 
                                  significance_level = significance_level)

    # relative observed benefit
    rel.obs.ben.mat[i, ] <- 
      observed_benefit_relative(status[group], w[group], significance_level = significance_level)
    
    # relative predicted benefit
    rel.pred.ben.mat[i, ] <- 
      predicted_benefit_inference(x = x, subset = group, 
                                  relative = TRUE, 
                                  benefits_risk = benefits_risk, 
                                  time_eval = time_eval, 
                                  significance_level = significance_level)
    
    # odds ratio
    or.mat[i, ] <- 
      odds_ratio(status[group], w[group], significance_level = significance_level)
    
  } # FOR
  
  return(list(predicted_benefit = list(absolute = abs.pred.ben.mat,
                                       relative = rel.pred.ben.mat),
              observed_benefit = list(absolute = abs.obs.ben.mat,
                                      relative = rel.obs.ben.mat),
              odds_ratio = or.mat,
              significance_level = significance_level,
              quantiles = colnames(quantile.groups),
              membership = quantile.groups))
} # FUN


#' Calculates group-level benefits from a GRF model, and the associated confidence intervals.
#' The returned benefits are the observed and predicted relative and absolute benefits as well as the odds ratio
#' 
#' @param x GRF model object
#' @param cutoffs the quantile cutoff points. Default is \code{c(0.25, 0.5, 0.75)}, which yields the quartiles.
#' @param baseline_risk The baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{x} is used.
#' @param significance_level the significance level. Default is 0.05.
#' 
#' @export
get_benefits_grf <- function(x, 
                             cutoffs = c(0.25, 0.5, 0.75),
                             baseline_risk = NULL,
                             significance_level = 0.05){
  
  # extract outcome and treatment status
  status <- x$inputs$status_bin
  w      <- x$inputs$w
  
  # specify baseline risk that shall be used for grouping
  if(is.null(baseline_risk)){
    baseline_risk <- as.numeric(x$risk$baseline)
  } # IF
  
  
  # group observations by their quantile of predicted baseline risk
  quantile.groups <- quantile_group_NoChecks(baseline_risk, cutoffs)
  
  ## calculate observed benefit and predicted benefit for each quantile group
  # initialize
  abs.obs.ben.mat           <- matrix(NA_real_, ncol(quantile.groups), 4L)
  colnames(abs.obs.ben.mat) <- c("estimate", "ci_lower", "ci_upper", "stderr")
  rownames(abs.obs.ben.mat) <- colnames(quantile.groups)
  abs.pred.ben.mat          <- abs.obs.ben.mat
  
  for(i in 1:ncol(quantile.groups)){
    
    group <- which(quantile.groups[,i])
    
    ## absolute observed benefit
    abs.obs.ben.mat[i,] <- 
      observed_benefit_absolute(status[group], w[group], significance_level = significance_level)
    
    ## absolute predicted benefit
    ate.group <- grf::average_treatment_effect(x$models, 
                                               subset = group)
    ate <- unname(ate.group["estimate"])
    se  <- unname(ate.group["std.err"]) 
    
    # quantile of the standard normal distribution
    z <- stats::qnorm(1 - significance_level/2)
    
    abs.pred.ben.mat[i, "estimate"] <- ate
    abs.pred.ben.mat[i, "stderr"]   <- se
    abs.pred.ben.mat[i, "ci_lower"] <- ate - z * se
    abs.pred.ben.mat[i, "ci_upper"] <- ate + z * se
    
  } # FOR group
  
  
  return(list(absolute_predicted_benefit = abs.pred.ben.mat, 
              absolute_observed_benefit = abs.obs.ben.mat, 
              significance_level = significance_level,
              quantiles = colnames(quantile.groups),
              membership = quantile.groups))
 
} # FUN





## helper function to get imputation-accounted absolute benefits
##
## @param x a 3-dimensional array as used in the function get.benefits_imputation.accounted()
## @param group.names the names of the quantile groups
## @param significance.level the significance level. Defaults to 5%
## @param relative. Logical. If TRUE, the quantity 'x' is a relative statistic. If FALSE, it is an absolute statistic.
##
## @return
benefits.imputed <- function(x, group.names, significance.level = 0.05, relative){
  
  # initialize
  m <- dim(x)[3]
  x[,"stderr",] <- x[,"stderr",]^2
  colnames(x)   <- c("estimate", "ci.lower", "ci.upper", "var")
  
  # get imputation-adjusted location and scale estimate
  T.hat <- rowMeans(x[,c("estimate"),])
  W.hat <- rowMeans(x[,c("var"),])
  B.hat <- rowSums((x[,c("estimate"),] - T.hat)^2) / (m-1)
  V.hat <- W.hat + (m+1) / m * B.hat
  
  # quantile of the standard normal distribution
  z <- stats::qnorm(1 - significance.level/2)
  
  # imputation-accounted standard deviation
  stderr <- sqrt(V.hat)
  
  # compute returned data frame depending on type of statistic
  if(relative){
    
    out <- data.frame(quantile = group.names, 
                      estimate = T.hat,
                      ci.lower = T.hat * exp(-z * stderr),
                      ci.upper = T.hat * exp(z * stderr),
                      stderr   = stderr)
    
  } else{
    
    out <- data.frame(quantile = group.names, 
                      estimate = T.hat,
                      ci.lower = T.hat - z * stderr,
                      ci.upper = T.hat + z * stderr,
                      stderr   = stderr)
    
  } # IF
  
  rownames(out) <- NULL
  return(out)
} # FUN



#' Calculates imputation-uncertainty-accounted group-level benefits from a prediction model, and the associated confidence intervals.
#' The returned benefits are the observed and predicted relative and absolute benefits as well as the odds ratio
#' 
#' @param pred.model.objs_imputed list of prediction model objects
#' @param cutoffs the quantile cutoff points. Default is c(0.25, 0.5, 0.75), which yields the quartiles.
#' @param risk.baseline A list of baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{pred.model.objs_imputed} is used.
#' @param significance.level the significance level. Default is 0.05.
#' 
#' @export
get.benefits_imputation.accounter <- function(pred.model.objs_imputed, 
                                              cutoffs = c(0.25, 0.5, 0.75),
                                              risk.baseline = NULL,
                                              significance.level = 0.05){
  
  # initialize 
  m                         <- length(pred.model.objs_imputed)
  group.names               <- intervals.quantile(cutoffs)
  rel.obs.ben.arr           <- array(NA_real_, dim = c(length(group.names), 4, m))
  colnames(rel.obs.ben.arr) <- c("estimate", "ci.lower", "ci.upper", "stderr")
  rownames(rel.obs.ben.arr) <- paste0("group.", 1:length(group.names))
  abs.obs.ben.arr           <- rel.obs.ben.arr
  abs.pred.ben.arr          <- rel.obs.ben.arr
  rel.pred.ben.arr          <- rel.obs.ben.arr
  or.arr                    <- rel.obs.ben.arr
  
  
  
  ## loop over the imputed datasets and account for the imputation uncertainty
  for(j in 1:m){
    
    # get j-th imputed model
    pred.model.obj <- pred.model.objs_imputed[[j]]
    
    # extract outcome and treatment status
    y <- pred.model.obj$inputs$y.prediction.timeframe
    w <- pred.model.obj$inputs$w
    
    # specify baseline risk that shall be used for grouping
    if(is.null(risk.baseline)){
      br <- pred.model.obj$risk$risk.baseline
    } else{
      br <- risk.baseline[[j]]
    }
    
    # group observations by their quantile of predicted baseline risk
    quantile.groups <- quantile_group(br, cutoffs)
    
    # get predicted benefit (relative and absolute)
    rel.pred.ben <- pred.model.obj$benefits$predicted.relative.benefit
    abs.pred.ben <- pred.model.obj$benefits$predicted.absolute.benefit
    
    ## calculate observed benefit and predicted benefit for each quantile group
    for(i in 1:ncol(quantile.groups)){
      group <- quantile.groups[,i]
      
      # absolute observed benefit
      abs.obs.ben.arr[i,,j] <- 
        absolute.observed.benefit(y[group], w[group], significance.level = significance.level)
      
      # absolute predicted benefit
      abs.pred.ben.arr[i,,j] <- 
        predicted.benefit.inference(abs.pred.ben[group], significance.level = significance.level)
      
      # relative observed benefit
      rel.obs.ben.arr[i,,j] <- 
        relative.observed.benefit(y[group], w[group], significance.level = significance.level)
      
      # relative predicted benefit
      rel.pred.ben.arr[i,,j] <- 
        predicted.benefit.inference(rel.pred.ben[group], significance.level = significance.level)
      
      # odds ratio
      or.arr[i,,j] <- 
        odds.ratio(y[group], w[group], significance.level = significance.level)
      
    } # FOR i 
    
  } # FOR imputed sets (j)
  
  
  # return imputation-accounted benefits
  return(list(absolute.observed.benefit  = benefits.imputed(x = abs.obs.ben.arr,
                                                            group.names = group.names, 
                                                            significance.level = significance.level,
                                                            relative = FALSE),
              absolute.predicted.benefit = benefits.imputed(x = abs.pred.ben.arr,
                                                            group.names = group.names, 
                                                            significance.level = significance.level,
                                                            relative = FALSE),
              relative.observed.benefit  = benefits.imputed(x = rel.obs.ben.arr,
                                                            group.names = group.names, 
                                                            significance.level = significance.level,
                                                            relative = TRUE),
              relative.predicted.benefit = benefits.imputed(x = rel.pred.ben.arr,
                                                            group.names = group.names, 
                                                            significance.level = significance.level,
                                                            relative = TRUE),
              odds.ratio                 = benefits.imputed(x = or.arr,
                                                            group.names = group.names, 
                                                            significance.level = significance.level,
                                                            relative = TRUE),
              significance.level         = significance.level,
              quantile.cutoff.points     = cutoffs))
} # FUN



#' Calculates imputation-uncertainty-accounted group-level benefits from a GRF model, and the associated confidence intervals.
#' The returned benefits are the observed and predicted absolute benefits as well as the odds ratio
#' 
#' @param grf.model.obj_imputed list of GRF model objects
#' @param cutoffs the quantile cutoff points. Default is c(0.25, 0.5, 0.75), which yields the quartiles.
#' @param risk.baseline A list of baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{grf.model.obj_imputed} is used.
#' @param significance.level the significance level. Default is 0.05.
#' 
#' @export 
get.benefits.grf_imputation.accounter <- function(grf.model.obj_imputed, 
                                                  cutoffs = c(0.25, 0.5, 0.75),
                                                  risk.baseline = NULL,
                                                  significance.level = 0.05){
  
  # initialize 
  m                         <- length(grf.model.obj_imputed)
  group.names               <- intervals.quantile(cutoffs)
  abs.obs.ben.arr           <- array(NA_real_, dim = c(length(group.names), 4, m))
  colnames(abs.obs.ben.arr) <- c("estimate", "ci.lower", "ci.upper", "stderr")
  rownames(abs.obs.ben.arr) <- paste0("group.", 1:length(group.names))
  abs.pred.ben.arr          <- abs.obs.ben.arr
  
  ## loop over the imputed datasets and account for the imputation uncertainty
  for(j in 1:m){
    
    # get j-th imputed model
    grf.model.obj <- grf.model.obj_imputed[[j]]
    
    # extract outcome and treatment status
    y <- grf.model.obj$inputs$y.prediction.timeframe
    w <- grf.model.obj$inputs$w
    
    # specify baseline risk that shall be used for grouping
    if(is.null(risk.baseline)){
      br <- grf.model.obj$risk$risk.baseline
    } else{
      br <- risk.baseline[[j]]
    }
    
    # group observations by their quantile of predicted baseline risk
    quantile.groups <- quantile_group(br, cutoffs)
    
    # get predicted benefit (absolute)
    abs.pred.ben <- grf.model.obj$benefits$predicted.absolute.benefit
    
    ## calculate observed benefit and predicted benefit for each quantile group
    for(i in 1:ncol(quantile.groups)){
      group <- quantile.groups[,i]
      
      # absolute observed benefit
      abs.obs.ben.arr[i,,j] <- 
        absolute.observed.benefit(y[group], w[group], significance.level = significance.level)
      
      # absolute predicted benefit
      abs.pred.ben.arr[i,,j] <- 
        predicted.benefit.inference(abs.pred.ben[group], significance.level = significance.level)
      
    } # FOR i 
    
  } # FOR imputed sets (j)
  
  
  # return imputation-accounted benefits
  return(list(absolute.observed.benefit  = benefits.imputed(x = abs.obs.ben.arr,
                                                            group.names = group.names, 
                                                            significance.level = significance.level,
                                                            relative = FALSE),
              absolute.predicted.benefit = benefits.imputed(x = abs.pred.ben.arr,
                                                            group.names = group.names, 
                                                            significance.level = significance.level,
                                                            relative = FALSE),
              significance.level         = significance.level,
              quantile.cutoff.points     = cutoffs))
} # FUN