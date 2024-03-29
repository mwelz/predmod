#' makes a calibration plot for a prediction model
#' 
#' @param x prediction model object, as returned by risk.modeling() or effect.modeling()
#' @param cutoffs the cutoff points of quantiles that shall be used for GATES grouping. Default is `c(0.25, 0.5, 0.75)`, which corresponds to the quartiles.
#' @param breaks Breaks along which to perform the grouping. If passed, overrules the grouping implied by \code{cutoffs}
#' @param relative logical. If `TRUE`, then relative benefits will be plotted. Default is `FALSE`
#' @param baseline_risk The baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{x} is used.
#' @param significance_level significance level for the confidence intervals. Default is 0.05
#' @param title optional title of the plot
#' @param xlim limits of x-axis
#' @param ylim limits of y-xcis
#' @param flip_sign logical. Shall the sign of the benefits be flipped?
#' @param newX Optional covariate matrix to calculate benefits with
#' @param newstatus Optional target variables to calculate benefits with
#' @param neww Optional treatment assignment variables to calculate benefits with
#' @param newz Optional linear predictors
#' @param shrunk TODO
#' 
#' @import ggplot2
#' 
#' @export
calibration_plot <- function( x,
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
  pb.means <- ob.means <- ob.means.ci.lo <- ob.means.ci.up <- NULL
  
  # get observed and predicted benefit by quantile group
  benefits <- get_benefits(x                  = x,
                           cutoffs            = cutoffs, 
                           breaks             = breaks,
                           baseline_risk      = baseline_risk,
                           odds_ratio         = FALSE,
                           significance_level = significance_level, 
                           newX               = newX,
                           newstatus          = newstatus,
                           neww               = neww,
                           newz               = newz, 
                           shrunk             = shrunk)
  
  # make sure risk quantile is in correct order
  risk.quantile <- factor(benefits$quantiles,
                          levels = benefits$quantiles)
  
  # adjust for relative and absolute benefit
  if(relative){
    
    if(!flip_sign)
    {
      df <- data.frame(pb.means = benefits$predicted_benefit$relative[,"estimate"],
                       ob.means = benefits$observed_benefit$relative[,"estimate"],
                       ob.means.ci.lo = benefits$observed_benefit$relative[,"ci_lower"],
                       ob.means.ci.up = benefits$observed_benefit$relative[,"ci_upper"],
                       risk.quantile = risk.quantile)
    } else
    {
      # flipping signs in a relative estimate means multiplying it by (-1) and adding one, 
      # that is, 1 - RR. We do the same for the CIs, which means that the originally 
      # "upper" bound of the CI becomes the lower bound and vice versa
      df <- data.frame(pb.means = 1.0 - benefits$predicted_benefit$relative[,"estimate"],
                       ob.means = 1.0 - benefits$observed_benefit$relative[,"estimate"],
                       ob.means.ci.lo = 1.0 - benefits$observed_benefit$relative[,"ci_upper"],
                       ob.means.ci.up = 1.0 - benefits$observed_benefit$relative[,"ci_lower"],
                       risk.quantile = risk.quantile)
    }
    
   
  } else{
    
    if(!flip_sign){
      
      df <- data.frame(pb.means = benefits$predicted_benefit$absolute[,"estimate"],
                       ob.means = benefits$observed_benefit$absolute[,"estimate"],
                       ob.means.ci.lo = benefits$observed_benefit$absolute[,"ci_lower"],
                       ob.means.ci.up = benefits$observed_benefit$absolute[,"ci_upper"],
                       risk.quantile = risk.quantile)
      
    } else{
      
      # when flipping signs, the originally "upper" bound of the CI becomes the lower bound
      # and vice versa
      df <- data.frame(pb.means = -benefits$predicted_benefit$absolute[,"estimate"],
                       ob.means = -benefits$observed_benefit$absolute[,"estimate"],
                       ob.means.ci.lo = -benefits$observed_benefit$absolute[,"ci_upper"],
                       ob.means.ci.up = -benefits$observed_benefit$absolute[,"ci_lower"],
                       risk.quantile = risk.quantile)
      
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


#' makes a calibration plot for a causal_forest model
#' 
#' @param x causal_forest model
#' @param cutoffs the cutoff points of quantiles that shall be used for GATES grouping. Default is `c(0.25, 0.5, 0.75)`, which corresponds to the quartiles.
#' @param breaks Breaks along which to perform the grouping. If passed, overrules the grouping implied by \code{cutoffs}
#' @param baseline_risk The baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{x} is used.
#' @param significance_level significance level for the confidence intervals. Default is 0.05
#' @param title optional title of the plot
#' @param xlim limits of x-axis
#' @param ylim limits of y-axis
#' @param flip_sign logical. Shall the sign of the benefits be flipped?
#' 
#' @import ggplot2
#' 
#' @export
calibration_plot_causal_forest <- function(x,
                                 cutoffs = c(0.25, 0.5, 0.75), 
                                 breaks = NULL,
                                 baseline_risk = NULL,
                                 significance_level = 0.05,
                                 title = NULL,
                                 xlim = NULL,
                                 ylim = NULL,
                                 flip_sign = FALSE){
  
  # appease the check (TODO: come up with better solution)
  pb.means <- ob.means <- ob.means.ci.lo <- ob.means.ci.up <- NULL
  
  # get observed and predicted benefit by quantile group
  benefits <- get_benefits_causal_forest(x = x, cutoffs = cutoffs, breaks = breaks,
                               baseline_risk = baseline_risk,
                               significance_level = significance_level)
  
  # make sure risk quantile is in correct order
  risk.quantile <- factor(benefits$quantiles,
                          levels = benefits$quantiles)
  
  if(!flip_sign){
    
    df <- data.frame(pb.means = benefits$absolute_predicted_benefit[,"estimate"],
                     ob.means = benefits$absolute_observed_benefit[,"estimate"],
                     ob.means.ci.lo = benefits$absolute_observed_benefit[,"ci_lower"],
                     ob.means.ci.up = benefits$absolute_observed_benefit[,"ci_upper"],
                     risk.quantile = risk.quantile)
    
  } else{
    
    df <- data.frame(pb.means = -benefits$absolute_predicted_benefit[,"estimate"],
                     ob.means = -benefits$absolute_observed_benefit[,"estimate"],
                     ob.means.ci.lo = -benefits$absolute_observed_benefit[,"ci_upper"],
                     ob.means.ci.up = -benefits$absolute_observed_benefit[,"ci_lower"],
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


group.static <- function(x, group.bounds){
  # helper function
  # if single value supplied: group is all obs with that particular value. 
  # if intervals are supplied: all obs within these bounds (closed set) are a group. Hence, intervals need to be non-overlapping. TODO: add option of having half-closed intervals. 
  # Note that rows that are all FALSE can be produced.
  
  group.mat <- matrix(NA, length(x), length(group.bounds))
  group.nam <- rep(NA_character_, length(group.bounds))
  
  for(i in 1:length(group.bounds)){
    bounds <- group.bounds[[i]]
    
    if(length(bounds) == 1){
      
      # case 1: a single value is supplied
      group.nam[i]   <- as.character(bounds)
      group.mat[, i] <- x == bounds
      
    } else if(length(bounds) == 2){
      
      # case 2: an interval is supplied
      group.nam[i]   <- paste0("[", bounds[1], ",", bounds[2], "]")
      group.mat[, i] <- (bounds[1] <= x) & (x <= bounds[2])
      
    } else stop("Incorrect input. Each element of 'group.bounds' needs to be a scalar or a 2-vector")
  }
  colnames(group.mat) <- group.nam
  return(group.mat)
} # FUN


order.intervals <- function(intervals, quantile.nam){
  # helper function to order intervals
  if(quantile.nam){
    is.first <- startsWith(intervals, "<")
    is.last  <- startsWith(intervals, ">")
    cut      <- as.numeric(gsub("[^0-9.-]", "",  gsub(",.*", "", intervals)))
    cut[is.first] <- -Inf
    cut[is.last]  <- Inf
    ord           <- order(cut, decreasing = FALSE)
  } else{
    ord <- order(as.numeric(substr(gsub(",.*", "", intervals), 2, 1e8)), 
                 decreasing = FALSE)
  }
  return(ord)
} # FUN


## make a subgroup plot of a variable x
## 
## @param pred.model.obj prediction model object, as returned by risk.modeling(), effect.modeling(), or grf.modeling()
## @param x a vector that shall be visualized
## @param relative TODO
## @param group.bounds TODO
## @param quantile.bounds TODO
## @param quantile.nam TODO
## @param risk.quantile.bounds TODO
## 
## @import ggplot2
## 
## @export
# subgroup.plot <- function(pred.model.obj, x, 
#                           relative = FALSE,
#                           group.bounds = NULL, 
#                           quantile.bounds = c(0.25, 0.5, 0.75), 
#                           quantile.nam = TRUE,
#                           risk.quantile.bounds = c(0.25, 0.5, 0.75)){
#   
#   # group.bounds <- list(-3.3, c(-3, -2), c(-2+0.0001, 2), c(2 + 0.0001, 3))
#   if(!is.null(group.bounds) & !is.null(quantile.bounds)) stop("Either group or quantile!")
#   if(!is.null(group.bounds) & !is.list(group.bounds)) stop("group.bounds needs to be a list")
#   
#   if(!is.null(group.bounds)){
#     x.group.mat <- group.static(x, group.bounds = group.bounds)
#   } else if(!is.null(quantile.bounds)){
#     x.group.mat <- quantile_group(x, cutoffs = quantile.bounds, names_quantile = quantile.nam)
#   } else stop("Please specify quantile bounds or group bounds")
#   
#   
#   # x-axis: risk group
#   risk.group.mat <- quantile_group(pred.model.obj$risk$risk.baseline, risk.quantile.bounds) 
#   risk.group <- rep(NA_character_, length(x))
#   for(nam in colnames(risk.group.mat)){
#     risk.group[which(risk.group.mat[,nam])] <- nam
#   }
#   
#   # drop the word 'quantile' and make sure x-axis is in logical order
#   risk.group <- factor(gsub( " .*$", "", risk.group))
#   lv <- levels(risk.group)
#   lv <- lv[order.intervals(lv, quantile.nam = TRUE)]
#   risk.group <- factor(risk.group, levels = lv)
#   
#   ## group by x's values
#   x.group <- rep(NA_character_, length(x))
#   for(nam in colnames(x.group.mat)){
#     x.group[which(x.group.mat[,nam])] <- nam
#   }
#   
#   if(quantile.nam){
#     # drop the word 'quantile' and make sure groups are in logical order
#     x.group <- factor(gsub( " .*$", "", x.group))
#     lv <- levels(x.group)
#     lv <- lv[order.intervals(lv, quantile.nam = TRUE)]
#     x.group <- factor(x.group, levels = lv)
#     leg.tit <- paste0("Quantile of ", deparse(substitute(x)))
#   } else{
#     x.group <- factor(gsub( " .*$", "", x.group))
#     lv <- levels(x.group)
#     lv <- lv[order.intervals(lv, quantile.nam = FALSE)]
#     x.group <- factor(x.group, levels = lv)
#     leg.tit <- paste0("Group of ", deparse(substitute(x)))
#   }
#   
#   # y-axis: predicted benefit
#   if(relative){
#     pred.ben <- pred.model.obj$benefits$predicted.relative.benefit
#   } else{
#     pred.ben <- pred.model.obj$benefits$predicted.absolute.benefit
#   }
#   
#   # prepare data frame
#   df <- stats::na.omit(data.frame(pred.ben, risk.group, x.group))
#   
#   ggplot(data = df, aes(x = risk.group, y = pred.ben, fill = x.group)) + 
#     geom_boxplot() +
#     theme_bw() +
#     xlab("Baseline risk quantile") +
#     ylab(ifelse(relative, "Predicted relative benefit", "Predicted absolute benefit")) +
#     scale_fill_discrete(name = leg.tit) +
#     theme(legend.position = "bottom")
# } # FUN
