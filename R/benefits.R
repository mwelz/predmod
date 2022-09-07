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
  
  # get absolute and relative risks
  get_predicted_benefits_NoChecks(risk_reg = risk_reg,
                                  risk_rev = risk_rev,
                                  w = w)
  
  
} # FUN


get_predicted_benefits_NoChecks <- function(risk_reg, risk_rev, w){
  
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
  
  y1 <- status[w==1]
  y0 <- status[w==0]
  
  ## account for scenarios where no t-tests can be performed: 
  # no variation or size less than two in either variable 
  # in this case, return zeros
  if(length(unique(y1)) < 2 |
     length(unique(y0)) < 2 |
     length(y1) < 2 | length(y0) < 2){
    
    out <- c(estimate = 0.0, 
             ci_lower = 0.0,
             ci_upper = 0.0,
             stderr   = 0.0)
  } else{
    
    ttest  <- stats::t.test(x = y1, y = y0, 
                            paired = FALSE, 
                            var.equal = FALSE) 
    aob    <- unname(ttest$estimate[1] - ttest$estimate[2])
    se.aob <- unname(ttest$stderr) 
    
    # quantile of the standard normal distribution
    z <- stats::qnorm(1 - significance_level/2)
    
    out <- c(estimate = aob, 
             ci_lower = aob - z * se.aob,
             ci_upper = aob + z * se.aob,
             stderr   = se.aob)
  } # IF
  return(out)
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
#' @param time_eval Time at which to evaluate the failure risk predictions.
#' @param significance_level The significance level, default is 0.05
#' 
#' @export
predicted_benefit_inference <- function(x, 
                                        subset = NULL, 
                                        relative = FALSE,
                                        time_eval = NULL,
                                        significance_level = 0.05){
  
  ate_obj <- average_treatment_effect_NoChecks(
                                      x = x, 
                                      subset = subset,
                                      relative = relative,
                                      time_eval = time_eval)
  
  se  <- unname(ate_obj["Std. Error"])
  ate <- unname(ate_obj["ATE"])

  # quantile of the standard normal distribution
  z <- stats::qnorm(1 - significance_level / 2)
  
  if(relative){
    cilo <- ate * exp(-z * se)
    ciup <- ate * exp( z * se)
  } else{
    cilo <- ate - z * se
    ciup <- ate + z * se
  } # IF
  
  # organize output to be consistent with other benefits functions
  return(c(estimate = ate, 
           ci_lower = cilo,
           ci_upper = ciup,
           stderr   = se))
} # FUN


#' Calculates group-level benefits from a prediction model, and the associated confidence intervals.
#' The returned benefits are the observed and predicted relative and absolute benefits as well as the odds ratio
#' 
#' @param x prediction model object
#' @param cutoffs the quantile cutoff points. Default is c(0.25, 0.5, 0.75), which yields the quartiles.
#' @param baseline_risk The baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{x} is used.
#' @param time_eval Time at which we evaluate the risk predictions.
#' @param odds_ratio Logical. If \code{TRUE}, odds ratios per quantile group will be computed. Default is \code{FALSE}.
#' @param significance_level the significance level. Default is 0.05.
#' @param status Optional target variables to calculate benefits with
#' @param w Optional treatment assignment variables to calculate benefits with
#' 
#' @export
get_benefits <- function(x, 
                         cutoffs = c(0.25, 0.5, 0.75),
                         baseline_risk = NULL,
                         time_eval = NULL,
                         odds_ratio = FALSE,
                         significance_level = 0.05, 
                         status = NULL,
                         w = NULL){
  
  status_null <- is.null(status)
  w_null <- is.null(w)
  baseline_null <- is.null(baseline_risk)

  if(is.null(status) && is.null(w))
  {
    
    # extract outcome and treatment status
    status <- x$inputs$status_bin
    w      <- x$inputs$w
    
  } else if(!is.null(status) && !is.null(w))
  {
    
    # nothing happens: take the supplied values as status and w, but run some tests
    stopifnot(is.numeric(status) && is.numeric(w))
    stopifnot(identical(length(status), length(w)))
    InputChecks_W(w)
    InputChecks_Y_binary(status)
    
  } else{
    stop(paste0("w and status must either be both NULL or both non-NULL"))
  } # IF
  
  ## if no baseline risk provided, take the baseline risk from x
  # but there might be no baseline risk in x if x didn't estimate a baseline risk.
  # in this case, throw an error.
  if(baseline_null)
  {
    br <- x$risk$baseline
    
    if(is.null(br)){
      
      # error if x$risk$baseline is NULL
      stop("Both the argument 'baseline_risk' and x$risk$baseline are NULL. ",
           "If x$risk$baseline = NULL, then baseline_risk cannot be NULL ",
           "as well. Please provide estimates of baseline risks via 'baseline_risk'.",
           call. = FALSE)
    } else{
      
      # x$risk$baseline contains baseline risk predictions. But they must be of same length as w and status!
      baseline_risk <- as.numeric(br)
    } # IF
  } # IF baseline_null
  
  ## we now have a baseline_risk object available, either implicitly obtained 
  # from x or explicitly as an argument. We now check for
  # equal length with w and status 
  if(!identical(length(baseline_risk), length(w)))
  {
    warning(paste0("You have not passed baseline_risk and ",
                   "the baseline_risk in x is of different length than ",
                   "w and status. This imbalance is no issue for the quantile ",
                   "grouping, but may be unintentional. So be careful here."))
  }

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
  if(odds_ratio){
    or.mat         <- rel.obs.ben.mat
  } else{
    or.mat         <- NULL
  } # IF
  
  
  for(i in 1:ncol(quantile.groups)){
    
    group <- which(quantile.groups[,i])

    # absolute observed benefit
    abs.obs.ben.mat[i, ] <- 
      observed_benefit_absolute(status = status[group], 
                                w = w[group], 
                                significance_level = significance_level)
    
    # absolute predicted benefit
    abs.pred.ben.mat[i, ] <- 
      predicted_benefit_inference(x = x, subset = group, 
                                  relative = FALSE, 
                                  time_eval = time_eval, 
                                  significance_level = significance_level)

    # relative observed benefit
    rel.obs.ben.mat[i, ] <- 
      observed_benefit_relative(status[group], w[group], significance_level = significance_level)
    
    # relative predicted benefit
    rel.pred.ben.mat[i, ] <- 
      predicted_benefit_inference(x = x, subset = group, 
                                  relative = TRUE, 
                                  time_eval = time_eval, 
                                  significance_level = significance_level)
    
    # odds ratio
    if(odds_ratio)
    {
      or.mat[i, ] <- 
        odds_ratio(status[group], w[group], significance_level = significance_level)
    } # IF
    
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
    if(class(x) == "grf_surv") stop("Baseline risk calculation in survival forests not yet implemented!")
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
