#' Calculates imputation-adjusted group-level benefits 
#' 
#' Calculates imputation-adjusted group-level benefits  from a prediction model, and the associated confidence intervals.The returned benefits are the observed and predicted relative and absolute benefits as well as the odds ratio.
#' @param x a list of prediction model object
#' @param cutoffs the quantile cutoff points. Default is c(0.25, 0.5, 0.75), which yields the quartiles.
#' @param breaks List of breaks along which to perform the grouping. If passed, overrules the grouping implied by \code{cutoffs}
#' @param baseline_risk A list of baseline risks that shall be used for grouping. If \code{NULL} (default), then the baseline risks in \code{x} are used.
#' @param significance_level the significance level. Default is 0.05.
#' @param newX Optional covariate matrix to calculate benefits with
#' @param newstatus Optional list of target variables to calculate benefits with
#' @param neww Optional list of treatment assignment variables to calculate benefits with
#' @param newz Optional linear predictors
#' @param shrunk TODO
#' @export
get_benefits_imputation <- function(x, 
                                    cutoffs = c(0.25, 0.5, 0.75),
                                    breaks = NULL,
                                    baseline_risk = NULL,
                                    significance_level = 0.05,
                                    newX = NULL,
                                    newstatus = NULL,
                                    neww = NULL,
                                    newz = NULL,
                                    shrunk = FALSE)
{
  # number of imputation runs
  m <- length(x)
  
  # sanity check
  if(!is.null(baseline_risk)) stopifnot(is.list(baseline_risk) & length(baseline_risk) == m)
  if(!is.null(newstatus)) stopifnot(is.list(newstatus) & length(newstatus) == m)
  if(!is.null(neww)) stopifnot(is.list(neww) & length(neww) == m)
  
  # initialize array of results
  arr <- array(data = NA_real_, 
               dim = c(length(cutoffs) + 1L, 2L, m),
               dimnames = list(NULL, c("estimate", "stderr"), NULL))
  
  # initialize list to store results in
  arr_ls <- list(pb_abs = arr,
                 pb_rel = arr,
                 ob_abs = arr,
                 ob_rel = arr) 
  #or     = arr) # don't calculate odds ration for now
  
  for(i in 1:m)
  {
    # get benefits
    ben <- get_benefits(x = x[[i]], 
                        cutoffs = cutoffs, 
                        breaks = breaks[[i]],
                        baseline_risk = baseline_risk[[i]], 
                        odds_ratio = FALSE, 
                        significance_level = significance_level, 
                        newX = newX[[i]], 
                        newstatus = newstatus[[i]], 
                        neww = neww[[i]], 
                        newz = newz[[i]], 
                        shrunk = shrunk)
    
    # assign results matrices
    arr_ls$pb_abs[,,i] <- ben$predicted_benefit$absolute[,c("estimate", "stderr")]
    arr_ls$pb_rel[,,i] <- ben$predicted_benefit$relative[,c("estimate", "stderr")]
    arr_ls$ob_abs[,,i] <- ben$observed_benefit$absolute[,c("estimate", "stderr")]
    arr_ls$ob_rel[,,i] <- ben$observed_benefit$relative[,c("estimate", "stderr")]
    # arr_ls$or[,,i]     <- ben$odds_ratio[,c("estimate", "stderr")]
    
  } # FOR
  
  # account for imputation uncertainty in absolute estimates
  pb_abs <- impaccount_regression_array(arr_ls$pb_abs,
                                        relative = FALSE, 
                                        significance_level = significance_level)
  ob_abs <- impaccount_regression_array(arr_ls$ob_abs,
                                        relative = FALSE, 
                                        significance_level = significance_level)
  
  # account for imputation uncertainty in relative estimates
  pb_rel <- impaccount_regression_array(arr_ls$pb_rel,
                                        relative = TRUE, 
                                        significance_level = significance_level)
  ob_rel <- impaccount_regression_array(arr_ls$ob_rel,
                                        relative = TRUE, 
                                        significance_level = significance_level)
  
  ## don't calculate odds ratio for now 
  # or     <- impaccount_regression_array(arr_ls$or,
  #                                       relative = TRUE, 
  #                                       significance_level = significance_level)
  
  
  
  ### adjust baseline risk 
  # take average over the baseline risk estimates
  if(is.null(baseline_risk))
  {
    br <- rowMeans(sapply(1:m, function(i) x[[i]]$risk$baseline))
  } else{
    br <- rowMeans(sapply(1:m, function(i) baseline_risk[[i]]))
  } # IF
  
  ## we now have a baseline_risk object available, either implicitly obtained 
  # from x or explicitly as an argument. We now check for
  # equal length with w and status 
  if(!is.null(neww) & !identical(length(br), length(neww[[1L]])))
  {
    warning(paste0("You have not passed baseline_risk and ",
                   "the baseline_risk in x is of different length than ",
                   "w and status. This imbalance is no issue for the quantile ",
                   "grouping, but may be unintentional. So be careful here."))
  }
  
  # group observations by their quantile of predicted baseline risk
  if(is.null(breaks))
  {
    quantile_groups <- quantile_group_NoChecks(br, cutoffs)
  } else{
    breaks0 <- rowMeans(sapply(1:m, function(i) breaks[[i]])) # average breaks value
    quantile_groups <- group_matrix(x = br, breaks = c(-Inf, breaks0, Inf))
  } # IF
  
  # group observations by their quantile of predicted baseline risk (imputation-adjusted)
  quantiles <- colnames(quantile_groups)
  
  # name the groups
  rownames(pb_abs) <- rownames(pb_rel) <- rownames(ob_abs) <- 
    rownames(ob_rel) <- quantiles
  
  # organize output
  out <- list(predicted_benefit = list(absolute = pb_abs,
                                       relative = pb_rel),
              observed_benefit = list(absolute = ob_abs,
                                      relative = ob_rel),
              #odds_ratio = or,
              significance_level = significance_level,
              risk_intervals = quantiles)
  return(out)
  
} # FUN



#' Calculates imputation-adjusted group-level benefits from a causal_forest model, and the associated confidence intervals.
#' The returned benefits are the observed and predicted relative and absolute benefits as well as the odds ratio
#' 
#' @param x A "causal_forest" object
#' @param cutoffs the quantile cutoff points. Default is \code{c(0.25, 0.5, 0.75)}, which yields the quartiles.
#' @param breaks List of breaks along which to perform the grouping. If passed, overrules the grouping implied by \code{cutoffs}
#' @param baseline_risk A list of baseline risks that shall be used for grouping. If \code{NULL} (default), then the baseline risks in \code{x} are used.
#' @param significance_level the significance level. Default is 0.05.
#' 
#' @export
get_benefits_causal_forest_imputation <- function(x, 
                                        cutoffs = c(0.25, 0.5, 0.75),
                                        breaks = NULL,
                                        baseline_risk = NULL,
                                        significance_level = 0.05)
{
  # number of imputation runs
  m <- length(x)
  
  # sanity check
  if(!is.null(baseline_risk)) stopifnot(is.list(baseline_risk) & length(baseline_risk) == m)
  
  # initialize array of results
  arr <- array(data = NA_real_, 
               dim = c(length(cutoffs) + 1L, 2L, m),
               dimnames = list(NULL, c("estimate", "stderr"), NULL))
  
  # initialize list to store results in
  arr_ls <- list(pb_abs = arr,
                 ob_abs = arr)
  
  for(i in 1:m)
  {
    # get benefits
    ben <- get_benefits_causal_forest(x = x[[i]], 
                            cutoffs = cutoffs, 
                            breaks = breaks[[i]],
                            baseline_risk = baseline_risk[[i]],
                            significance_level = significance_level)
    
    # assign results matrices
    arr_ls$pb_abs[,,i] <- ben$absolute_predicted_benefit[,c("estimate", "stderr")]
    arr_ls$ob_abs[,,i] <- ben$absolute_observed_benefit[,c("estimate", "stderr")]
    
  } # FOR
  
  # account for imputation uncertainty in absolute estimates
  pb_abs <- impaccount_regression_array(arr_ls$pb_abs,
                                        relative = FALSE, 
                                        significance_level = significance_level)
  ob_abs <- impaccount_regression_array(arr_ls$ob_abs,
                                        relative = FALSE, 
                                        significance_level = significance_level)
  
  ### adjust baseline risk 
  # take average over the baseline risk estimates
  if(is.null(baseline_risk))
  {
    br <- rowMeans(sapply(1:m, function(i) x[[i]]$risk$baseline))
  } else{
    br <- rowMeans(sapply(1:m, function(i) baseline_risk[[i]]))
  } # IF
  
  ## we now have a baseline_risk object available, either implicitly obtained 
  # from x or explicitly as an argument. We now check for
  # equal length with w and status 
  if(!identical(length(br), length(x[[1]]$inputs$w)))
  {
    warning(paste0("You have not passed baseline_risk and ",
                   "the baseline_risk in x is of different length than ",
                   "w and status. This imbalance is no issue for the quantile ",
                   "grouping, but may be unintentional. So be careful here."))
  }
  
  # group observations by their quantile of predicted baseline risk (imputation-adjusted)
  if(is.null(breaks))
  {
    quantile_groups <- quantile_group_NoChecks(br, cutoffs)
  } else{
    breaks0 <- rowMeans(sapply(1:m, function(i) breaks[[i]])) # average breaks value
    quantile_groups <- group_matrix(x = br, breaks = c(-Inf, breaks0, Inf))
  } # IF
  
  quantiles <- colnames(quantile_groups)
  
  # name the groups
  rownames(pb_abs) <- rownames(ob_abs) <- quantiles
  
  # organize output
  out <- list(absolute_predicted_benefit = pb_abs,
              absolute_observed_benefit = ob_abs,
              significance_level = significance_level,
              risk_intervals = quantiles)
  return(out)
  
} # FUN
