
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
}

#' makes a calibration plot for a prediction model
#' 
#' @param x prediction model object, as returned by risk.modeling() or effect.modeling()
#' @param cutoffs the cutoff points of quantiles that shall be used for GATES grouping. Default is `c(0.25, 0.5, 0.75)`, which corresponds to the quartiles.
#' @param relative logical. If `TRUE`, then relative benefits will be plotted. Default is `FALSE`
#' @param baseline_risk The baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{pred.model.obj} is used.
#' @param benefits_risk Logical. If \code{TRUE}, then the risk-based benefits are used (only applicable to survival models). Default is \code{FALSE}.
#' @param time_eval Only applicable if \code{benefits_risk = TRUE}. Time at which we evaluate the risk predictions.
#' @param significance_level significance level for the confidence intervals. Default is 0.05
#' @param title optional title of the plot
#' @param xlim limits of x-axis
#' @param ylim limits of y-xcis
#' @param flip_sign logical. Shall the sign of the benefits be flipped?
#' 
#' @import ggplot2
#' 
#' @export
calibration_plot <- function( x,
                              cutoffs = c(0.25, 0.5, 0.75), 
                              relative = FALSE,
                              baseline_risk = NULL,
                              benefits_risk = FALSE,
                              time_eval = NULL,
                              significance_level = 0.05,
                              title = NULL,
                              xlim = NULL,
                              ylim = NULL,
                              flip_sign = FALSE){
  
  # appease the check (TODO: come up with better solution)
  pb.means <- ob.means <- ob.means.ci.lo <- ob.means.ci.up <- NULL
  
  # get observed and predicted benefit by quantile group
  benefits <- get_benefits(x                  = x,
                           cutoffs            = cutoffs,
                           baseline_risk      = baseline_risk,
                           benefits_risk      = benefits_risk,
                           time_eval          = time_eval,
                           significance_level = significance_level)
  
  # make sure risk quantile is in correct order
  risk.quantile <- factor(benefits$quantiles,
                          levels = benefits$quantiles)
  
  # adjust for relative and absolute benefit
  if(relative){
    
    df <- data.frame(pb.means = benefits$predicted_benefit$relative[,"estimate"],
                     ob.means = benefits$observed_benefit$relative[,"estimate"],
                     ob.means.ci.lo = benefits$observed_benefit$relative[,"ci_lower"],
                     ob.means.ci.up = benefits$observed_benefit$relative[,"ci_upper"],
                     risk.quantile = risk.quantile)
  } else{
    
    if(!flip_sign){
      
      df <- data.frame(pb.means = benefits$predicted_benefit$absolute[,"estimate"],
                       ob.means = benefits$observed_benefit$absolute[,"estimate"],
                       ob.means.ci.lo = benefits$observed_benefit$absolute[,"ci_lower"],
                       ob.means.ci.up = benefits$observed_benefit$absolute[,"ci_upper"],
                       risk.quantile = risk.quantile)
      
    } else{
      
      df <- data.frame(pb.means = -benefits$predicted_benefit$absolute[,"estimate"],
                       ob.means = -benefits$observed_benefit$absolute[,"estimate"],
                       ob.means.ci.lo = -benefits$observed_benefit$absolute[,"ci_lower"],
                       ob.means.ci.up = -benefits$observed_benefit$absolute[,"ci_upper"],
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


#' makes a calibration plot for a GRF model
#' 
#' @param grf.model.obj GRF model oobject as returned by grf.modeling()
#' @param cutoffs the cutoff points of quantiles that shall be used for GATES grouping. Default is `c(0.25, 0.5, 0.75)`, which corresponds to the quartiles.
#' @param risk.baseline The baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{pred.model.obj} is used.
#' @param significance.level significance level for the confidence intervals. Default is 0.05
#' @param title optional title of the plot
#' @param xlim limits of x-axis
#' @param ylim limits of y-xcis
#' @param flip.sign.of.absolute.benefit logical. Shall the sign of the benefits be flipped?
#' 
#' @import ggplot2
#' 
#' @export
calibration.plot.grf <- function(grf.model.obj,
                                 cutoffs = c(0.25, 0.5, 0.75), 
                                 risk.baseline = NULL,
                                 significance.level = 0.05,
                                 title = NULL,
                                 xlim = NULL,
                                 ylim = NULL,
                                 flip.sign.of.absolute.benefit = FALSE){
  
  # appease the check (TODO: come up with better solution)
  pb.means <- ob.means <- ob.means.ci.lo <- ob.means.ci.up <- NULL
  
  # get observed and predicted benefit by quantile group
  benefits <- get.benefits.grf(grf.model.obj, 
                               cutoffs = cutoffs,
                               risk.baseline = risk.baseline,
                               significance.level = significance.level)
  
  # make sure risk quantile is in correct order
  risk.quantile <- factor(benefits$absolute.observed.benefit$quantile,
                          levels = benefits$absolute.observed.benefit$quantile)

  # retrieve data
  if(!flip.sign.of.absolute.benefit){
    
    df <- data.frame(pb.means = benefits$absolute.predicted.benefit$estimate,
                     ob.means = benefits$absolute.observed.benefit$estimate,
                     ob.means.ci.lo = benefits$absolute.observed.benefit$ci.lower,
                     ob.means.ci.up = benefits$absolute.observed.benefit$ci.upper,
                     risk.quantile = risk.quantile)
    
  } else{
    
    df <- data.frame(pb.means = -benefits$absolute.predicted.benefit$estimate,
                     ob.means = -benefits$absolute.observed.benefit$estimate,
                     ob.means.ci.lo = -benefits$absolute.observed.benefit$ci.lower,
                     ob.means.ci.up = -benefits$absolute.observed.benefit$ci.upper,
                     risk.quantile = risk.quantile)
    
  }
  
  
  if(is.null(title)) title <- "Calibration plot of absolute benefit"
  
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



#' make a subgroup plot of a variable x
#' 
#' @param pred.model.obj prediction model object, as returned by risk.modeling(), effect.modeling(), or grf.modeling()
#' @param x a vector that shall be visualized
#' @param relative TODO
#' @param group.bounds TODO
#' @param quantile.bounds TODO
#' @param quantile.nam TODO
#' @param risk.quantile.bounds TODO
#' 
#' @import ggplot2
#' 
#' @export
subgroup.plot <- function(pred.model.obj, x, 
                          relative = FALSE,
                          group.bounds = NULL, 
                          quantile.bounds = c(0.25, 0.5, 0.75), 
                          quantile.nam = TRUE,
                          risk.quantile.bounds = c(0.25, 0.5, 0.75)){
  
  # group.bounds <- list(-3.3, c(-3, -2), c(-2+0.0001, 2), c(2 + 0.0001, 3))
  if(!is.null(group.bounds) & !is.null(quantile.bounds)) stop("Either group or quantile!")
  if(!is.null(group.bounds) & !is.list(group.bounds)) stop("group.bounds needs to be a list")
  
  if(!is.null(group.bounds)){
    x.group.mat <- group.static(x, group.bounds = group.bounds)
  } else if(!is.null(quantile.bounds)){
    x.group.mat <- quantile_group(x, cutoffs = quantile.bounds, names_quantile = quantile.nam)
  } else stop("Please specify quantile bounds or group bounds")
  
  
  # x-axis: risk group
  risk.group.mat <- quantile_group(pred.model.obj$risk$risk.baseline, risk.quantile.bounds) 
  risk.group <- rep(NA_character_, length(x))
  for(nam in colnames(risk.group.mat)){
    risk.group[which(risk.group.mat[,nam])] <- nam
  }
  
  # drop the word 'quantile' and make sure x-axis is in logical order
  risk.group <- factor(gsub( " .*$", "", risk.group))
  lv <- levels(risk.group)
  lv <- lv[order.intervals(lv, quantile.nam = TRUE)]
  risk.group <- factor(risk.group, levels = lv)
  
  ## group by x's values
  x.group <- rep(NA_character_, length(x))
  for(nam in colnames(x.group.mat)){
    x.group[which(x.group.mat[,nam])] <- nam
  }
  
  if(quantile.nam){
    # drop the word 'quantile' and make sure groups are in logical order
    x.group <- factor(gsub( " .*$", "", x.group))
    lv <- levels(x.group)
    lv <- lv[order.intervals(lv, quantile.nam = TRUE)]
    x.group <- factor(x.group, levels = lv)
    leg.tit <- paste0("Quantile of ", deparse(substitute(x)))
  } else{
    x.group <- factor(gsub( " .*$", "", x.group))
    lv <- levels(x.group)
    lv <- lv[order.intervals(lv, quantile.nam = FALSE)]
    x.group <- factor(x.group, levels = lv)
    leg.tit <- paste0("Group of ", deparse(substitute(x)))
  }
  
  # y-axis: predicted benefit
  if(relative){
    pred.ben <- pred.model.obj$benefits$predicted.relative.benefit
  } else{
    pred.ben <- pred.model.obj$benefits$predicted.absolute.benefit
  }
  
  # prepare data frame
  df <- stats::na.omit(data.frame(pred.ben, risk.group, x.group))
  
  ggplot(data = df, aes(x = risk.group, y = pred.ben, fill = x.group)) + 
    geom_boxplot() +
    theme_bw() +
    xlab("Baseline risk quantile") +
    ylab(ifelse(relative, "Predicted relative benefit", "Predicted absolute benefit")) +
    scale_fill_discrete(name = leg.tit) +
    theme(legend.position = "bottom")
} # FUN


#' makes a calibration plot for a prediction model while accounting for imputation uncertainty
#' 
#' @param pred.model.objs_imputed prediction model object, as returned by risk.modeling() or effect.modeling() TODO: use name of imputation functions
#' @param cutoffs the cutoff points of quantiles that shall be used for GATES grouping. Default is c(0.25, 0.5, 0.75), which corresponds to the quartiles.
#' @param relative logical. If TRUE, then relative benefits will be plotted. Default is FALSE
#' @param risk.baseline A list of baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{pred.model.objs_imputed} is used.
#' @param significance.level significance level for the confidence intervals. Default is 0.05
#' @param title optional title of the plot
#' @param xlim limits of x-axis
#' @param ylim limits of y-xcis
#' @param flip.sign.of.absolute.benefit logical. Shall the sign of the benefits be flipped?
#' 
#' @export
calibration.plot_imputation.accounter <- function(pred.model.objs_imputed,
                                                  cutoffs = c(0.25, 0.5, 0.75), 
                                                  relative = FALSE,
                                                  risk.baseline = NULL,
                                                  significance.level = 0.05,
                                                  title = NULL,
                                                  xlim = NULL,
                                                  ylim = NULL,
                                                  flip.sign.of.absolute.benefit = FALSE){
  
  # appease the check (TODO: come up with better solution)
  pb.means <- ob.means <- ob.means.ci.lo <- ob.means.ci.up <- NULL
  
  # get observed and predicted benefit by quantile group
  benefits <- get.benefits_imputation.accounter(pred.model.objs_imputed, 
                                                cutoffs = cutoffs,
                                                risk.baseline = risk.baseline,
                                                significance.level = significance.level)
  
  # get average risk quantile 
  if(is.null(risk.baseline)){
    
    rb <- rowMeans(sapply(1:length(pred.model.objs_imputed), function(i){
      
      pred.model.objs_imputed[[i]]$risk$risk.baseline
      
    }))
    
  } else{
    
    rb <- rowMeans(sapply(1:length(pred.model.objs_imputed), function(i){
      
      risk.baseline[[i]]
      
    }))
    
  } # IF
  
  qg <- quantile_group(rb, cutoffs = cutoffs)
  risk.quantile <- factor(colnames(qg), levels = colnames(qg))
  
  # adjust for relative and absolute benefit
  if(relative){
    
    df <- data.frame(pb.means = benefits$relative.predicted.benefit$estimate,
                     ob.means = benefits$relative.observed.benefit$estimate,
                     ob.means.ci.lo = benefits$relative.observed.benefit$ci.lower,
                     ob.means.ci.up = benefits$relative.observed.benefit$ci.upper,
                     risk.quantile = risk.quantile)
  } else{
    
    if(!flip.sign.of.absolute.benefit){
      
      df <- data.frame(pb.means = benefits$absolute.predicted.benefit$estimate,
                       ob.means = benefits$absolute.observed.benefit$estimate,
                       ob.means.ci.lo = benefits$absolute.observed.benefit$ci.lower,
                       ob.means.ci.up = benefits$absolute.observed.benefit$ci.upper,
                       risk.quantile = risk.quantile)
      
    } else{
      
      df <- data.frame(pb.means = -benefits$absolute.predicted.benefit$estimate,
                       ob.means = -benefits$absolute.observed.benefit$estimate,
                       ob.means.ci.lo = -benefits$absolute.observed.benefit$ci.lower,
                       ob.means.ci.up = -benefits$absolute.observed.benefit$ci.upper,
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


#' makes a calibration plot for a GRF model while accounting for imputation uncertainty
#' 
#' @param grf.model.obj_imputed GRF model object, as returned by risk.modeling() or effect.modeling() TODO: use name of imputation functions
#' @param cutoffs the cutoff points of quantiles that shall be used for GATES grouping. Default is `c(0.25, 0.5, 0.75)`, which corresponds to the quartiles.
#' @param risk.baseline A list of baseline risk that shall be used for grouping. If \code{NULL} (default), then the baseline risk in \code{pred.model.objs_imputed} is used.
#' @param relative logical. If `TRUE`, then relative benefits will be plotted. Default is `FALSE`
#' @param significance.level significance level for the confidence intervals. Default is 0.05
#' @param title optional title of the plot
#' @param xlim limits of x-axis
#' @param ylim limits of y-xcis
#' @param flip.sign.of.absolute.benefit logical. Shall the sign of the benefits be flipped?
#' 
#' @export
calibration.plot.grf_imputation.accounter <- function(grf.model.obj_imputed,
                                                      cutoffs = c(0.25, 0.5, 0.75), 
                                                      risk.baseline = NULL,
                                                      relative = FALSE,
                                                      significance.level = 0.05,
                                                      title = NULL,
                                                      xlim = NULL,
                                                      ylim = NULL,
                                                      flip.sign.of.absolute.benefit = FALSE){
  
  # appease the check (TODO: come up with better solution)
  pb.means <- ob.means <- ob.means.ci.lo <- ob.means.ci.up <- NULL
  
  # get observed and predicted benefit by quantile group
  benefits <- get.benefits.grf_imputation.accounter(grf.model.obj_imputed, 
                                                    cutoffs = cutoffs, 
                                                    risk.baseline = risk.baseline,
                                                    significance.level = significance.level)
  
  # get average risk quantile 
  if(is.null(risk.baseline)){
    
    rb <- rowMeans(sapply(1:length(grf.model.obj_imputed), function(i){
      
      grf.model.obj_imputed[[i]]$risk$risk.baseline
      
    }))
    
  } else{
    
    rb <- rowMeans(sapply(1:length(grf.model.obj_imputed), function(i){
      
      risk.baseline[[i]]
      
    }))
    
  } # IF
  
  qg <- quantile_group(rb, cutoffs = cutoffs)
  risk.quantile <- factor(colnames(qg), levels = colnames(qg))
  
  # adjust for relative and absolute benefit
  if(!flip.sign.of.absolute.benefit){
    
    df <- data.frame(pb.means = benefits$absolute.predicted.benefit$estimate,
                     ob.means = benefits$absolute.observed.benefit$estimate,
                     ob.means.ci.lo = benefits$absolute.observed.benefit$ci.lower,
                     ob.means.ci.up = benefits$absolute.observed.benefit$ci.upper,
                     risk.quantile = risk.quantile)
    
  } else{
    
    df <- data.frame(pb.means = -benefits$absolute.predicted.benefit$estimate,
                     ob.means = -benefits$absolute.observed.benefit$estimate,
                     ob.means.ci.lo = -benefits$absolute.observed.benefit$ci.lower,
                     ob.means.ci.up = -benefits$absolute.observed.benefit$ci.upper,
                     risk.quantile = risk.quantile)
    
  } # IF
  

  if(is.null(title)) title <- "Calibration plot of absolute benefit"
  
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