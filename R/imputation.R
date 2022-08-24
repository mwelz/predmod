# x is a 3D-array of regression outputs on the 3rd dimension
impaccount_regression_array <- function(x, 
                                        relative = FALSE, 
                                        significance_level = 0.05){
  
  # since x carries 'stderr' in 2nd dimension, but we need the variance, we need to adjust
  x[,2L,] <- x[,2L,]^2L
  
  # get the relevant statistics
  m     <- dim(x)[3L]
  ave   <- apply(x, c(1L,2L), mean)
  T_hat <- ave[, 1L] # location estimate
  W_hat <- ave[, 2L] # within-variance estimate
  
  # between-variance estimate
  B_hat <- rowSums(sapply(1:m, 
                          function(i) (x[,1L,i] - T_hat)^2))/(m-1)
  
  # pooled variance estimate
  V_hat <- W_hat + (m+1)/m * B_hat
  
  # get confidence interval
  stderr <- sqrt(V_hat)
  z      <- stats::qnorm(1 - significance_level / 2)
  
  if(relative){
    cilo <- T_hat * exp(-z * stderr)
    ciup <- T_hat * exp(z * stderr)
  } else{
    cilo <- T_hat - z * stderr
    ciup <- T_hat + z * stderr
  } # IF
  
  # organize output to be consistent with other benefits functions
  out <- cbind(estimate = T_hat,
               ci_lower = cilo,
               ci_upper = ciup,
               stderr = stderr)
  return(out)
  
} # FUN


#' Calculates imputation-adjusted group-level benefits from a prediction model, and the associated confidence intervals.
#' The returned benefits are the observed and predicted relative and absolute benefits as well as the odds ratio
#' 
#' @param x a list of prediction model object
#' @param cutoffs the quantile cutoff points. Default is c(0.25, 0.5, 0.75), which yields the quartiles.
#' @param baseline_risk A list of baseline risks that shall be used for grouping. If \code{NULL} (default), then the baseline risks in \code{x} are used.
#' @param time_eval Time at which we evaluate the risk predictions.
#' @param significance_level the significance level. Default is 0.05.
#' 
#' @export
get_benefits_imputation <- function(x, 
                                    cutoffs = c(0.25, 0.5, 0.75),
                                    baseline_risk = NULL,
                                    time_eval = NULL,
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
                 pb_rel = arr,
                 ob_abs = arr,
                 ob_rel = arr) 
                 #or     = arr) # don't calculate odds ration for now
  
  for(i in 1:m)
  {
    # get benefits
    ben <- get_benefits(x = x[[i]], 
                        cutoffs = cutoffs,
                        baseline_risk = baseline_risk[[i]], 
                        time_eval = time_eval, 
                        odds_ratio = FALSE,
                        significance_level = significance_level)
    
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
  
  # group observations by their quantile of predicted baseline risk (imputation-adjusted)
  quantile_groups <- quantile_group_NoChecks(br, cutoffs)
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


#' Calculates imputation-adjusted group-level benefits from a GRF model, and the associated confidence intervals.
#' The returned benefits are the observed and predicted relative and absolute benefits as well as the odds ratio
#' 
#' @param x GRF model object
#' @param cutoffs the quantile cutoff points. Default is \code{c(0.25, 0.5, 0.75)}, which yields the quartiles.
#' @param baseline_risk A list of baseline risks that shall be used for grouping. If \code{NULL} (default), then the baseline risks in \code{x} are used.
#' @param significance_level the significance level. Default is 0.05.
#' 
#' @export
get_benefits_grf_imputation <- function(x, 
                                        cutoffs = c(0.25, 0.5, 0.75),
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
    ben <- get_benefits_grf(x = x[[i]], 
                            cutoffs = cutoffs,
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

  # organize output
  out <- list(absolute_predicted_benefit = pb_abs,
              absolute_observed_benefit = ob_abs,
              significance_level = significance_level)
  return(out)
  
} # FUN


#' makes an imputation-accounted calibration plot for imputed prediction models
#' 
#' @param x A list of prediction model objects
#' @param cutoffs the cutoff points of quantiles that shall be used for GATES grouping. Default is `c(0.25, 0.5, 0.75)`, which corresponds to the quartiles.
#' @param relative logical. If `TRUE`, then relative benefits will be plotted. Default is `FALSE`
#' @param baseline_risk A list of baseline risks that shall be used for grouping. If \code{NULL} (default), then the baseline risks as in \code{x} are used.
#' @param time_eval Time at which we evaluate the risk predictions.
#' @param significance_level significance level for the confidence intervals. Default is 0.05
#' @param title optional title of the plot
#' @param xlim limits of x-axis
#' @param ylim limits of y-xcis
#' @param flip_sign logical. Shall the sign of the benefits be flipped?
#' 
#' @import ggplot2
#' 
#' @export
calibration_plot_imputation <- function(x,
                                        cutoffs = c(0.25, 0.5, 0.75), 
                                        relative = FALSE,
                                        baseline_risk = NULL,
                                        time_eval = NULL,
                                        significance_level = 0.05,
                                        title = NULL,
                                        xlim = NULL,
                                        ylim = NULL,
                                        flip_sign = FALSE){
  
  # appease the check (TODO: come up with better solution)
  pb.means <- ob.means <- ob.means.ci.lo <- ob.means.ci.up <- NULL
  
  # get observed and predicted benefit by quantile group (imputation-adjusted)
  benefits <- get_benefits_imputation(x                  = x,
                                      cutoffs            = cutoffs,
                                      baseline_risk      = baseline_risk,
                                      time_eval          = time_eval,
                                      significance_level = significance_level)
  
  
  
  # make sure risk quantile is in correct order
  risk_quantile <- factor(benefits$risk_intervals,
                          levels = benefits$risk_intervals)
  
  # adjust for relative and absolute benefit
  if(relative){
    
    df <- data.frame(pb.means = benefits$predicted_benefit$relative[,"estimate"],
                     ob.means = benefits$observed_benefit$relative[,"estimate"],
                     ob.means.ci.lo = benefits$observed_benefit$relative[,"ci_lower"],
                     ob.means.ci.up = benefits$observed_benefit$relative[,"ci_upper"],
                     risk.quantile = risk_quantile)
  } else{
    
    if(!flip_sign){
      
      df <- data.frame(pb.means = benefits$predicted_benefit$absolute[,"estimate"],
                       ob.means = benefits$observed_benefit$absolute[,"estimate"],
                       ob.means.ci.lo = benefits$observed_benefit$absolute[,"ci_lower"],
                       ob.means.ci.up = benefits$observed_benefit$absolute[,"ci_upper"],
                       risk.quantile = risk_quantile)
      
    } else{
      
      df <- data.frame(pb.means = -benefits$predicted_benefit$absolute[,"estimate"],
                       ob.means = -benefits$observed_benefit$absolute[,"estimate"],
                       ob.means.ci.lo = -benefits$observed_benefit$absolute[,"ci_lower"],
                       ob.means.ci.up = -benefits$observed_benefit$absolute[,"ci_upper"],
                       risk.quantile = risk_quantile)
      
    }
  } # IF
  
  
  if(is.null(title)){
    title <- paste0("Calibration plot of ", 
                    ifelse(relative, "relative ", "absolute "), "benefit")
  } # IF
  
  ggplot(mapping = aes(x = pb.means,
                       y = ob.means, color = risk_quantile), data = df) +
    geom_point() +
    geom_errorbar(mapping = aes(ymin = ob.means.ci.lo,
                                ymax = ob.means.ci.up)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = ifelse(relative, "Predicted relative benefit", "Predicted absolute benefit"),
         y = ifelse(relative, "Observed relative benefit", "Observed absolute benefit")) +
    theme_light() +
    ggtitle(title) +
    theme(legend.position = "bottom") 
  
} # FUN




#' makes a calibration plot for a GRF model and adjusts for imputation uncertainty
#' 
#' @param x A list of grf models
#' @param cutoffs the cutoff points of quantiles that shall be used for GATES grouping. Default is `c(0.25, 0.5, 0.75)`, which corresponds to the quartiles.
#' @param baseline_risk A list of baseline risks that shall be used for grouping. If \code{NULL} (default), then the baseline risks in \code{x} are used.
#' @param significance_level significance level for the confidence intervals. Default is 0.05
#' @param title optional title of the plot
#' @param xlim limits of x-axis
#' @param ylim limits of y-axis
#' @param flip_sign logical. Shall the sign of the benefits be flipped?
#' 
#' @import ggplot2
#' 
#' @export
calibration_plot_grf_imputation <- function(x,
                                            cutoffs = c(0.25, 0.5, 0.75), 
                                            baseline_risk = NULL,
                                            significance_level = 0.05,
                                            title = NULL,
                                            xlim = NULL,
                                            ylim = NULL,
                                            flip_sign = FALSE){
  
  # appease the check (TODO: come up with better solution)
  pb.means <- ob.means <- ob.means.ci.lo <- ob.means.ci.up <- NULL
  
  # get observed and predicted benefit by quantile group (imputation-adjusted)
  benefits <- get_benefits_grf_imputation(x                  = x,
                                          cutoffs            = cutoffs,
                                          baseline_risk      = baseline_risk,
                                          significance_level = significance_level)
  
  # get imputation-adjusted baseline risk
  m <- length(x)
  if(is.null(baseline_risk))
  {
    br <- rowMeans(sapply(1:m, function(i) x[[i]]$risk$baseline))
  } else{
    br <- rowMeans(sapply(1:m, function(i) baseline_risk[[i]]))
  } # IF
  
  # group observations by their quantile of predicted baseline risk (imputation-adjusted)
  quantile.groups <- quantile_group_NoChecks(br, cutoffs)
  quantiles <- colnames(quantile.groups)
  
  # make sure risk quantile is in correct order
  risk.quantile <- factor(quantiles,
                          levels = quantiles)
  
  # flip sign if requested
  if(!flip_sign){
    
    df <- data.frame(pb.means = benefits$absolute_predicted_benefit[,"estimate"],
                     ob.means = benefits$absolute_observed_benefit[,"estimate"],
                     ob.means.ci.lo = benefits$absolute_observed_benefit[,"ci_lower"],
                     ob.means.ci.up = benefits$absolute_observed_benefit[,"ci_upper"],
                     risk.quantile = risk.quantile)
    
  } else{
    
    df <- data.frame(pb.means = -benefits$absolute_predicted_benefit[,"estimate"],
                     ob.means = -benefits$absolute_observed_benefit[,"estimate"],
                     ob.means.ci.lo = -benefits$absolute_observed_benefit[,"ci_lower"],
                     ob.means.ci.up = -benefits$absolute_observed_benefit[,"ci_upper"],
                     risk.quantile = risk.quantile)
    
  } # IF
  
  if(is.null(title)){
    title <- "Calibration plot of absolute benefit"
  } # IF
  
  ggplot(mapping = aes(x = pb.means,
                       y = ob.means, color = risk.quantile), data = df) +
    geom_point() +
    geom_errorbar(mapping = aes(ymin = ob.means.ci.lo,
                                ymax = ob.means.ci.up)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = "Predicted absolute benefit",
         y = "Observed absolute benefit") +
    theme_light() +
    ggtitle(title) +
    theme(legend.position = "bottom") 
  
} # FUN


# C is list of concordance estimates, SE a list of associated standard errors
impaccount_concordance <- function(C, SE)
{
  m <- length(C)
  T_hat <- mean(sapply(1:m, function(i) C[[i]]))
  W_hat <- mean(sapply(1:m, function(i) SE[[i]]^2)) # square since it's a variance
  B_hat <- sum(sapply(1:m, function(i) (C[[i]] - T_hat)^2 )) / (m - 1)
  V_hat <- W_hat + (m+1)/m * B_hat
  
  return(list(estimate = T_hat, stderr = sqrt(V_hat)))
  
} # FUN


#' Account for imputation uncertainty in risk models
#' 
#' @param x a list of effect model objects (imputed)
#' @return A list of imputation-accounted stuff
#' 
#' @export
impaccount_risk_model <- function(x)
{
  # number of imputation runs
  m <- length(x)
  
  ## prepare lists
  # benfits
  benefits <- list(absolute = NULL, relative = NULL)
  benefits$absolute <- rowMeans(sapply(1:m, function(i) x[[i]]$benefits$absolute))
  benefits$relative <- rowMeans(sapply(1:m, function(i) x[[i]]$benefits$relative))
  
  # coefficients
  coefficients <- list(baseline = NULL, stage2 = NULL)
  
  # account for the case that baseline risk is NULL
  if(all(sapply(1:m, function(i) !is.null(x[[i]]$coefficients$baseline ))))
  {
    cf <- Matrix::Matrix(rowMeans(sapply(1:m, 
                                         function(i) as.matrix(x[[i]]$coefficients$baseline))),
                         sparse = TRUE)
    dimnames(cf) <- dimnames(x[[1]]$coefficients$baseline)
    coefficients$baseline <- cf
  }
  
  decisions <- c("accepted", "rejected")
  
  for(decision in decisions)
  {
    p <- nrow(x[[1]]$coefficients$stage2[[decision]])
    arr <- array(NA_real_, dim = c(p, 4L, m))
    for(i in 1:m) arr[,,i] <- x[[i]]$coefficients$stage2[[decision]]
    tmp <- impaccount_regression_array(x = arr, relative = FALSE)
    z <- tmp[, "estimate"] / tmp[, "stderr"]
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    cf <- cbind(tmp[, "estimate"], tmp[, "stderr"], z, p)
    dimnames(cf) <- dimnames(x[[1]]$coefficients$stage2[[decision]])
    coefficients$stage2[[decision]] <- cf
  } # FOR decision
  
  
  # risk
  risk <- list(baseline = NULL, regular = NULL, counterfactual = NULL)
  risk$regular        <- as.matrix(rowMeans(sapply(1:m, function(i) x[[i]]$risk$regular)))
  risk$counterfactual <- as.matrix(rowMeans(sapply(1:m, function(i) x[[i]]$risk$counterfactual)))
  
  if(all(sapply(1:m, function(i) !is.null(x[[i]]$risk$baseline ))))
  {
    risk$baseline     <- as.matrix(rowMeans(sapply(1:m, function(i) x[[i]]$risk$baseline)))
  }
  
  
  ## concordance
  # concordance <- list(outcome_baseline = NULL,
  #                     outcome = NULL,
  #                     benefit = NULL)
  # 
  # c_logi  <- all(sapply(1:m, function(i) !is.null(x[[i]]$concordance$outcome_baseline$estimate )))
  # se_logi <- all(sapply(1:m, function(i) !is.null(x[[i]]$concordance$outcome_baseline$stderr )))
  # if(c_logi & se_logi){
  #   concordance$outcome_baseline <- 
  #     impaccount_concordance(C  = lapply(1:m, function(i) x[[i]]$concordance$outcome_baseline$estimate),
  #                            SE = lapply(1:m, function(i) x[[i]]$concordance$outcome_baseline$stderr))
  # }
  # 
  # concordance$outcome <- 
  #   impaccount_concordance(C  = lapply(1:m, function(i) x[[i]]$concordance$outcome$estimate),
  #                          SE = lapply(1:m, function(i) x[[i]]$concordance$outcome$stderr))
  # concordance$benefit <- 
  #   impaccount_concordance(C  = lapply(1:m, function(i) x[[i]]$concordance$benefit$estimate),
  #                          SE = lapply(1:m, function(i) x[[i]]$concordance$benefit$stderr))
  
  
  return(list(benefits = benefits,
              coefficients = coefficients,
              risk = risk))
} # FOR


#' Account for imputation uncertainty in effect models
#' 
#' @param x a list of effect model objects (imputed)
#' @return A list of imputation-accounted stuff
#' 
#' @export
impaccount_effect_model <- function(x)
{
  # number of imputation runs
  m <- length(x)
  
  ## prepare lists
  # benfits
  benefits <- list(absolute = NULL, relative = NULL)
  benefits$absolute <- rowMeans(sapply(1:m, function(i) x[[i]]$benefits$absolute))
  benefits$relative <- rowMeans(sapply(1:m, function(i) x[[i]]$benefits$relative))
  
  # coefficients
  coefficients <- list(baseline = NULL, full = NULL, reduced = NULL)
  
  # account for the case that baseline risk is NULL
  if(all(sapply(1:m, function(i) !is.null(x[[i]]$coefficients$baseline ))))
  {
    cf <- Matrix::Matrix(rowMeans(sapply(1:m, 
                                         function(i) as.matrix(x[[i]]$coefficients$baseline))),
                         sparse = TRUE)
    dimnames(cf) <- dimnames(x[[1]]$coefficients$baseline)
    coefficients$baseline <- cf
  }
  
  coefficients$full <- 
    Matrix::Matrix(rowMeans(sapply(1:m, 
                                   function(i) as.matrix(x[[i]]$coefficients$full))),
                   sparse = TRUE)
  
  # jointly retained variables
  var_names <- Reduce(intersect, lapply(1:m, function(i) rownames(x[[i]]$coefficients$reduced)))
  
  arr <- array(NA_real_, dim = c(length(var_names), 4L, m))
  for(i in 1:m) arr[,,i] <- x[[i]]$coefficients$reduced[var_names,]
  tmp <- impaccount_regression_array(x = arr, relative = FALSE)
  z <- tmp[, "estimate"] / tmp[, "stderr"]
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  cf <- cbind(tmp[, "estimate"], tmp[, "stderr"], z, p)
  colnames(cf) <- colnames(x[[1]]$coefficients$reduced)
  rownames(cf) <- var_names
  coefficients$reduced <- cf
  
  # risk
  risk <- list(baseline = NULL, regular = NULL, counterfactual = NULL)
  risk$regular        <- as.matrix(rowMeans(sapply(1:m, function(i) x[[i]]$risk$regular)))
  risk$counterfactual <- as.matrix(rowMeans(sapply(1:m, function(i) x[[i]]$risk$counterfactual)))
  
  if(all(sapply(1:m, function(i) !is.null(x[[i]]$risk$baseline ))))
  {
    risk$baseline     <- as.matrix(rowMeans(sapply(1:m, function(i) x[[i]]$risk$baseline)))
  }
  
  # concordance
  concordance <- list(outcome_baseline = NULL,
                      outcome = NULL,
                      benefit = NULL)
  
  concordance$outcome <- 
    impaccount_concordance(C  = lapply(1:m, function(i) x[[i]]$concordance$outcome$estimate),
                           SE = lapply(1:m, function(i) x[[i]]$concordance$outcome$stderr))
  concordance$benefit <- 
    impaccount_concordance(C  = lapply(1:m, function(i) x[[i]]$concordance$benefit$estimate),
                           SE = lapply(1:m, function(i) x[[i]]$concordance$benefit$stderr))
  
  c_logi  <- all(sapply(1:m, function(i) !is.null(x[[i]]$concordance$outcome_baseline$estimate )))
  se_logi <- all(sapply(1:m, function(i) !is.null(x[[i]]$concordance$outcome_baseline$stderr )))
  if(c_logi & se_logi){
    concordance$outcome_baseline <- 
      impaccount_concordance(C  = lapply(1:m, function(i) x[[i]]$concordance$outcome_baseline$estimate),
                             SE = lapply(1:m, function(i) x[[i]]$concordance$outcome_baseline$stderr))
  }
  
  
  return(list(benefits = benefits,
              coefficients = coefficients,
              risk = risk, 
              concordance = concordance))
} # FOR
