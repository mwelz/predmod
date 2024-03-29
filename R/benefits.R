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
#' @param significance_level The significance level, default is 0.05
#' @param neww Optional treatment assignment variables to calculate benefits with
#' @param newX Optional covariate matrix to calculate benefits with
#' @param newz Optional linear predictors
#' @param shrunk TODO
#' 
#' @export
predicted_benefit_inference <- function(x, 
                                        subset = NULL, 
                                        relative = FALSE,
                                        significance_level = 0.05,
                                        neww = NULL, 
                                        newX = NULL,
                                        newz = NULL,
                                        shrunk = FALSE){
  
  ate_obj <- average_treatment_effect(
                                      x = x, 
                                      subset = subset,
                                      relative = relative,
                                      neww = neww, 
                                      newX = newX,
                                      newz = newz,
                                      shrunk = shrunk)
  
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
#' @param breaks Breaks along which to perform the grouping. If passed, overrules the grouping implied by \code{cutoffs}
#' @param baseline_risk The baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{x} is used.
#' @param odds_ratio Logical. If \code{TRUE}, odds ratios per quantile group will be computed. Default is \code{FALSE}.
#' @param significance_level the significance level. Default is 0.05.
#' @param newX Optional covariate matrix to calculate benefits with
#' @param newstatus Optional target variables to calculate benefits with
#' @param neww Optional treatment assignment variables to calculate benefits with
#' @param newz Optional linear predictors
#' @param shrunk TODO
#' 
#' @export
get_benefits <- function(x,
                         cutoffs = c(0.25, 0.5, 0.75),
                         breaks = NULL,
                         baseline_risk = NULL,
                         odds_ratio = FALSE,
                         significance_level = 0.05, 
                         newX = NULL,
                         newstatus = NULL,
                         neww = NULL,
                         newz = NULL,
                         shrunk = FALSE){
  
  stopifnot(inherits(x = x, what = c("risk_model_crss", "effect_model_crss", "causal_forest")))
  if(!is.null(cutoffs) && !is.null(breaks))
  {
    message("Both cutoffs and breaks were passed. Breaks will be used for grouping")
  }
  
  ## the input checks of newX, neww, newz will be performed in the ATE functions,
  ## so no need to do them here. Instead, check for consistency with baseline risk
  NULLw <- is.null(neww)
  NULLy <- is.null(newstatus)
  
  if(NULLw && NULLy)
  {
    w0      <- x$inputs$w
    status0 <- x$inputs$status_bin
  } else if(!NULLw && !NULLy)
  {
    w0      <- neww
    status0 <- newstatus

  } else
  {
    stop("neww and newstatus must be either both NULL or both non-NULL")
  }
  
  # input check and return adjusted baseline risks
  br <- InputChecks_get_benefits_and_return_br(x = x, 
                                               status = status0, 
                                               w = w0, 
                                               baseline_risk = baseline_risk)
  
  # group observations by their quantile of predicted baseline risk
  if(is.null(breaks))
  {
    quantile.groups <- quantile_group_NoChecks(br, cutoffs)
  } else{
    quantile.groups <- group_matrix(x = br, breaks = c(-Inf, breaks, Inf))
  } # IF
  
  
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
  
  
  for(i in seq_len(ncol(quantile.groups))){
    
    group <- which(quantile.groups[,i])

    # absolute observed benefit 
    abs.obs.ben.mat[i, ] <- 
      observed_benefit_absolute(status = status0[group], 
                                w = w0[group], 
                                significance_level = significance_level) # w cannot be NULL here
    
    # absolute predicted benefit
    abs.pred.ben.mat[i, ] <- 
      predicted_benefit_inference(x = x, subset = group, 
                                  relative = FALSE, 
                                  significance_level = significance_level, 
                                  neww = neww, newX = newX, newz = newz, shrunk = shrunk) # keep them as passed

    # relative observed benefit
    rel.obs.ben.mat[i, ] <- 
      observed_benefit_relative(status = status0[group], 
                                w = w0[group], 
                                significance_level = significance_level)  # w cannot be NULL here
    
    # relative predicted benefit
    rel.pred.ben.mat[i, ] <- 
      predicted_benefit_inference(x = x, subset = group, 
                                  relative = TRUE, 
                                  significance_level = significance_level,
                                  neww = neww, newX = newX, newz = newz, shrunk = shrunk) # keep them as passed
    
    # odds ratio
    if(odds_ratio)
    {
      or.mat[i, ] <- 
        odds_ratio(status0[group], w0[group], significance_level = significance_level)
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


#' Calculates group-level benefits from a causal_forest model, and the associated confidence intervals.
#' The returned benefits are the observed and predicted relative and absolute benefits as well as the odds ratio
#' 
#' @param x causal_forest model object
#' @param cutoffs the quantile cutoff points. Default is \code{c(0.25, 0.5, 0.75)}, which yields the quartiles.
#' @param breaks Breaks along which to perform the grouping. If passed, overrules the grouping implied by \code{cutoffs}
#' @param baseline_risk The baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{x} is used.
#' @param significance_level the significance level. Default is 0.05.
#' 
#' @export
get_benefits_causal_forest <- function(x, 
                             cutoffs = c(0.25, 0.5, 0.75),
                             breaks = NULL,
                             baseline_risk = NULL,
                             significance_level = 0.05){
  
  stopifnot(inherits(x, what = "causal_forest"))
  if(!is.null(cutoffs) && !is.null(breaks))
  {
    message("Both cutoffs and breaks were passed. Breaks will be used for grouping")
  }
  
  # extract outcome and treatment status
  status   <- x$inputs$status_bin
  w        <- x$inputs$w
  predbens <- x$benefits$absolute
  
  # specify baseline risk that shall be used for grouping
  # input check and return adjusted baseline risks
  baseline_risk <- InputChecks_get_benefits_and_return_br(
    x = x, 
    status = x$inputs$status_bin, 
    w = x$inputs$w, 
    baseline_risk = baseline_risk)
  
  # group observations by their quantile of predicted baseline risk
  if(is.null(breaks))
  {
    quantile.groups <- quantile_group_NoChecks(baseline_risk, cutoffs)
  } else{
    quantile.groups <- group_matrix(x = baseline_risk, breaks = c(-Inf, breaks, Inf))
  } # IF
  
  ## calculate observed benefit and predicted benefit for each quantile group
  # initialize
  abs.obs.ben.mat           <- matrix(NA_real_, ncol(quantile.groups), 4L)
  colnames(abs.obs.ben.mat) <- c("estimate", "ci_lower", "ci_upper", "stderr")
  rownames(abs.obs.ben.mat) <- colnames(quantile.groups)
  abs.pred.ben.mat          <- abs.obs.ben.mat
  
  # quantile of the standard normal distribution
  z <- stats::qnorm(1 - significance_level/2)
  
  for(i in 1:ncol(quantile.groups)){
    
    group <- which(quantile.groups[,i])
    
    ## absolute observed benefit
    abs.obs.ben.mat[i,] <- 
      observed_benefit_absolute(status[group], w[group], significance_level = significance_level)
    
    ## absolute predicted benefit. Don't do the AIPW-ATE of the GRF package, as this may mess up the order of the predbens
    # think of the AIPW-ATE as a weighted mean of predbens: if the weights are nonequal, the sorted (weighted) predbens may
    # give rise to a different size ordering than the unweighted predbens
    t <- stats::t.test(x = predbens[group])
    ate <- unname(t$estimate)
    se  <- t$stderr

    # ate.group <- grf::average_treatment_effect(x$models, 
    #                                            subset = group)
    # ate <- unname(ate.group["estimate"])
    # se  <- unname(ate.group["std.err"]) 
    
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


# helper function for get_beenfits(). x is any predmod object
# returns baseline risk vector
InputChecks_get_benefits_and_return_br <- 
  function(x, 
           status, # non-NULL
           w, # non-NULL
           baseline_risk = NULL)
{
 
  # nothing happens: take the supplied values as status and w, but run some tests
  stopifnot(is.numeric(status) && is.numeric(w))
  stopifnot(identical(length(status), length(w)))
  InputChecks_W(w)
  InputChecks_Y_binary(status)
  
  ## if no baseline risk provided, take the baseline risk from x
  # but there might be no baseline risk in x if x didn't estimate a baseline risk.
  # in this case, throw an error.
  baseline_null <- is.null(baseline_risk)
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
      
      # x$risk$baseline contains baseline risk predictions
      baseline_risk <- as.numeric(br)
    } # IF
  } # IF baseline_null
  
  ## we now have a baseline_risk object available, either implicitly obtained 
  # from x or explicitly as an argument. We now check for
  # equal length with w and status 
  if(!identical(length(baseline_risk), length(w)))
  {
    stop(paste0("You have not passed baseline_risk and ",
                   "the baseline_risk in x is of different length than ",
                   "the length of treatment assignment and status that you",
                   " wish to use for the calculation of benefits"))
  } # IF
  
  return(baseline_risk)
  
} # FUN
