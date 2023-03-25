#' makes an imputation-accounted calibration plot for imputed prediction models
#' 
#' @param x A list of prediction model objects
#' @param cutoffs the cutoff points of quantiles that shall be used for GATES grouping. Default is `c(0.25, 0.5, 0.75)`, which corresponds to the quartiles.
#' @param breaks List of breaks along which to perform the grouping. If passed, overrules the grouping implied by \code{cutoffs}
#' @param relative logical. If `TRUE`, then relative benefits will be plotted. Default is `FALSE`
#' @param baseline_risk A list of baseline risks that shall be used for grouping. If \code{NULL} (default), then the baseline risks as in \code{x} are used.
#' @param significance_level significance level for the confidence intervals. Default is 0.05
#' @param title optional title of the plot
#' @param xlim limits of x-axis
#' @param ylim limits of y-xcis
#' @param flip_sign logical. Shall the sign of the benefits be flipped?
#' @param newstatus Optional list of target variables to calculate benefits with
#' @param neww Optional list of treatment assignment variables to calculate benefits with
#' 
#' @import ggplot2
#' 
#' @export
calibration_plot_imputation <- function(x,
                                        cutoffs = c(0.25, 0.5, 0.75), 
                                        breaks = NULL,
                                        relative = FALSE,
                                        baseline_risk = NULL,
                                        significance_level = 0.05,
                                        title = NULL,
                                        xlim = NULL,
                                        ylim = NULL,
                                        flip_sign = FALSE, 
                                        newX = NULL,
                                        newstatus = NULL,
                                        neww = NULL,
                                        newz = NULL,
                                        shrunk = FALSE){
  
  # appease the check (TODO: come up with better solution)
  pb.means <- ob.means <- ob.means.ci.lo <- ob.means.ci.up <- risk.quantile <- NULL
  
  # get observed and predicted benefit by quantile group (imputation-adjusted)
  benefits <- get_benefits_imputation(x                  = x,
                                      cutoffs            = cutoffs,
                                      breaks             = breaks,
                                      baseline_risk      = baseline_risk,
                                      significance_level = significance_level,
                                      newX               = newX, 
                                      newstatus          = newstatus, 
                                      neww               = neww, 
                                      newz               = newz,
                                      shrunk             = shrunk)
  
  
  
  # make sure risk quantile is in correct order
  risk_quantile <- factor(benefits$risk_intervals,
                          levels = benefits$risk_intervals)
  
  # adjust for relative and absolute benefit
  if(relative){
    
    if(!flip_sign)
    {
      df <- data.frame(pb.means = benefits$predicted_benefit$relative[,"estimate"],
                       ob.means = benefits$observed_benefit$relative[,"estimate"],
                       ob.means.ci.lo = benefits$observed_benefit$relative[,"ci_lower"],
                       ob.means.ci.up = benefits$observed_benefit$relative[,"ci_upper"],
                       risk.quantile = risk_quantile)
    } else
    {
      # flipping signs in a relative estimate means multiplying it by (-1) and adding one, 
      # that is, 1 - RR. We do the same for the CIs, which means that the originally 
      # "upper" bound of the CI becomes the lower bound and vice versa
      df <- data.frame(pb.means = 1.0 - benefits$predicted_benefit$relative[,"estimate"],
                       ob.means = 1.0 - benefits$observed_benefit$relative[,"estimate"],
                       ob.means.ci.lo = 1.0 - benefits$observed_benefit$relative[,"ci_upper"],
                       ob.means.ci.up = 1.0 - benefits$observed_benefit$relative[,"ci_lower"],
                       risk.quantile = risk_quantile)
    }
    
    
  } else{
    
    if(!flip_sign){
      
      df <- data.frame(pb.means = benefits$predicted_benefit$absolute[,"estimate"],
                       ob.means = benefits$observed_benefit$absolute[,"estimate"],
                       ob.means.ci.lo = benefits$observed_benefit$absolute[,"ci_lower"],
                       ob.means.ci.up = benefits$observed_benefit$absolute[,"ci_upper"],
                       risk.quantile = risk_quantile)
      
    } else{
      
      # when flipping signs, the originally "upper" bound of the CI becomes the lower bound
      # and vice versa
      df <- data.frame(pb.means = -benefits$predicted_benefit$absolute[,"estimate"],
                       ob.means = -benefits$observed_benefit$absolute[,"estimate"],
                       ob.means.ci.lo = -benefits$observed_benefit$absolute[,"ci_upper"],
                       ob.means.ci.up = -benefits$observed_benefit$absolute[,"ci_lower"],
                       risk.quantile = risk_quantile)
      
    }
  } # IF
  
  
  if(is.null(title)){
    title <- paste0("Calibration plot of ", 
                    ifelse(relative, "relative ", "absolute "), "benefit")
  } # IF
  
  ggplot(mapping = aes(x = pb.means,
                       y = ob.means, color = risk.quantile), data = df) +
    geom_point() +
    geom_errorbar(mapping = aes(ymin = ob.means.ci.lo,
                                ymax = ob.means.ci.up)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_vline(xintercept = 0, linetype = 3) +
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
#' @param breaks List of breaks along which to perform the grouping. If passed, overrules the grouping implied by \code{cutoffs}
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
                                            breaks = NULL,
                                            baseline_risk = NULL,
                                            significance_level = 0.05,
                                            title = NULL,
                                            xlim = NULL,
                                            ylim = NULL,
                                            flip_sign = FALSE){
  
  # appease the check (TODO: come up with better solution)
  pb.means <- ob.means <- ob.means.ci.lo <- ob.means.ci.up <- risk.quantile <- NULL
  
  # get observed and predicted benefit by quantile group (imputation-adjusted)
  benefits <- get_benefits_grf_imputation(x                  = x,
                                          cutoffs            = cutoffs,
                                          breaks             = breaks,
                                          baseline_risk      = baseline_risk,
                                          significance_level = significance_level)
  
  # make sure risk quantile is in correct order
  risk_quantile <- factor(benefits$risk_intervals,
                          levels = benefits$risk_intervals)
  
  # flip sign if requested
  if(!flip_sign){
    
    df <- data.frame(pb.means = benefits$absolute_predicted_benefit[,"estimate"],
                     ob.means = benefits$absolute_observed_benefit[,"estimate"],
                     ob.means.ci.lo = benefits$absolute_observed_benefit[,"ci_lower"],
                     ob.means.ci.up = benefits$absolute_observed_benefit[,"ci_upper"],
                     risk.quantile = risk_quantile)
    
  } else{
    
    df <- data.frame(pb.means = -benefits$absolute_predicted_benefit[,"estimate"],
                     ob.means = -benefits$absolute_observed_benefit[,"estimate"],
                     ob.means.ci.lo = -benefits$absolute_observed_benefit[,"ci_upper"],
                     ob.means.ci.up = -benefits$absolute_observed_benefit[,"ci_lower"],
                     risk.quantile = risk_quantile)
    
  } # IF
  
  if(is.null(title)){
    title <- "Calibration plot of absolute benefit"
  } # IF
  
  ggplot(mapping = aes(x = pb.means,
                       y = ob.means, color = risk.quantile), data = df) +
    geom_point() +
    geom_errorbar(mapping = aes(ymin = ob.means.ci.lo,
                                ymax = ob.means.ci.up)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_vline(xintercept = 0, linetype = 3) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = "Predicted absolute benefit",
         y = "Observed absolute benefit") +
    theme_light() +
    ggtitle(title) +
    theme(legend.position = "bottom") 
  
} # FUN

