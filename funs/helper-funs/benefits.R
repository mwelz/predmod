#' groups the observations according to their quantiles
#' 
#' @param x the numeric vector to be grouped
#' @param cutoffs the quantile cutoff points. Default is c(0.25, 0.5, 0.75), which yields the quartiles.
#' @param quantile.nam shall the names of the quantiles or their cutoff values be used for the group keys? Default is TRUE (quantile names)
#' 
#' @export
quantile.group <- function(x, cutoffs = c(0.25, 0.5, 0.75), quantile.nam = TRUE){
  # cutoffs are the quantile cutoffs (like c(0.25, 0.5, 0.75))
  q         <- quantile(x, cutoffs)
  q         <- c(-Inf, q, Inf)
  groups    <- as.character(cut(x, breaks = q, include.lowest = FALSE, right = TRUE, dig.lab = 3))
  group.nam <- unique(groups)
  group.nam <- group.nam[order(
    as.numeric(substr(sub("\\,.*", "", group.nam), 2, stop = 1e8L)), 
    decreasing = FALSE)] # ensure the order is correct
  group.mat <- matrix(NA, length(x), length(group.nam))
  nam       <- rep(NA, length(group.nam))
  
  for(j in 1:length(group.nam)){
    if(j == 1){
      nam[j] <- paste0("<=", 100*cutoffs[j], "% quantile")
    } else if (j == length(group.nam)){
      nam[j] <- paste0(">", 100*cutoffs[j-1], "% quantile")
    } else{
      nam[j] <- paste0("(", 100*cutoffs[j-1], ",", 100*cutoffs[j], "]% quantile")
    }
    group.mat[,j] <- groups == group.nam[j]
  }
  
  if(quantile.nam){
    colnames(group.mat) <- nam
  } else{
    colnames(group.mat) <- group.nam
  }
  return(group.mat)
}


#' calculates the absolute observed benefit, which corresponds to the difference in mean(y[W=w]), and the associated confidence interval.
#' @param y vector of binary outcomes
#' @param w vector of binary treatment status
#' @param significance.level the significance level. Default is 0.05.
#'
#' @export
absolute.observed.benefit <- function(y, w, significance.level = 0.05){
  
  ttest  <- t.test(y[w==1], y[w==0], paired = FALSE, var.equal = FALSE) 
  aob    <- unname(ttest$estimate[1] - ttest$estimate[2])
  se.aob <- unname(ttest$stderr) 
  
  # quantile of the standard normal distribution
  z <- qnorm(1 - significance.level/2)
  
  return(c(estimate = aob, 
           ci.lower = aob - z * se.aob,
           ci.upper = aob + z * se.aob,
           stderr   = se.aob))
}


#' calculates the relative observed benefit, which corresponds to a risk ratio, and the associated confidence interval
#' @param y vector of binary outcomes
#' @param w vector of binary treatment status
#' @param significance.level the significance level. Default is 0.05.
#'
#' @export
relative.observed.benefit <- function(y, w, significance.level = 0.05){
  
  #' contingency table with the following structure
  #'      Y=1 Y=0
  #'  W=1  a   b
  #'  W=0  c   d
  a <- sum(w == 1 & y == 1)
  b <- sum(w == 1 & y == 0)
  c <- sum(w == 0 & y == 1)
  d <- sum(w == 0 & y == 0)
  
  # estimates for the probabilities
  p1.hat <- a / (a + b) # estimates Pr(Y = 1 | W = 1)
  p0.hat <- c / (c + d) # estimates Pr(Y = 1 | W = 0)
  
  # empirical risk ratio (the relative observed benefit) and its standard error
  rob    <- p1.hat / p0.hat
  se.rob <- sqrt(b / (a * (a + b)) + d / (c * (c + d)))
  
  # quantile of the standard normal distribution
  z <- qnorm(1 - significance.level/2)
  
  return(c(estimate = rob, 
           ci.lower = rob * exp(-z * se.rob),
           ci.upper = rob * exp(z * se.rob),
           stderr   = se.rob))
}


#' calculates the odds ratio, and the associated confidence interval. 
#' @param w vector of binary treatment status
#' @param significance.level the significance level. Default is 0.05.
#'
#' @export
odds.ratio <- function(y, w, significance.level = 0.05){

  #' contingency table with the following structure
  #'      Y=1 Y=0
  #'  W=1  a   b
  #'  W=0  c   d
  a <- sum(w == 1 & y == 1)
  b <- sum(w == 1 & y == 0)
  c <- sum(w == 0 & y == 1)
  d <- sum(w == 0 & y == 0)
  
  # estimates for the probabilities
  p1.hat <- a / (a + b) # estimates Pr(Y = 1 | W = 1)
  p0.hat <- c / (c + d) # estimates Pr(Y = 1 | W = 0)
  
  # empirical odds ratio and its standard error
  or <- (p1.hat / (1 - p1.hat) ) / (p0.hat / (1 - p0.hat))
  se.or <- sqrt(1/a + 1/b + 1/c + 1/d)
  
  # quantile of the standard normal distribution
  z <- qnorm(1 - significance.level/2)
  
  return(c(estimate = or, 
           ci.lower = or * exp(-z * se.or),
           ci.upper = or * exp(z * se.or),
           stderr   = se.or))
}


#' performs inference on the per-individual predicted benefits
#' 
#' @param predicted.benefits a vector of per-individual predicted benefits (either absolute of relative)
#' @param significance.level the significance level. Default is 0.05.
#' 
#' @export
predicted.benefit.inference <- function(predicted.benefits, significance.level = 0.05){
  
  ttest <- t.test(predicted.benefits)
  pb    <- unname(ttest$estimate)
  se.pb <- unname(ttest$stderr)
  
  # quantile of the standard normal distribution
  z <- qnorm(1 - significance.level/2)
  
  return(c(estimate = pb, 
           ci.lower = pb - z * se.pb,
           ci.upper = pb + z * se.pb,
           stderr   = se.pb))
}


#' Calculates group-level benefits from a prediction model, and the associated confidence intervals.
#' The returned benefits are the observed and predicted relative and absolute benefits as well as the odds ratio
#' 
#' @param pred.model.object prediction model object
#' @param cutoffs the quantile cutoff points. Default is c(0.25, 0.5, 0.75), which yields the quartiles.
#' @param significance.level the significance level. Default is 0.05.
#' 
#' @export
get.benefits <- function(pred.model.obj, 
                         cutoffs = c(0.25, 0.5, 0.75),
                         significance.level = 0.05){
  
  # extract outcome and treatment status
  y <- pred.model.obj$inputs$y.prediction.timeframe
  w <- pred.model.obj$inputs$w
  
  # group observations by their quantile of predicted baseline risk
  quantile.groups <- quantile.group(pred.model.obj$risk$risk.baseline, cutoffs)
  
  # get predicted benefit (relative and absolute)
  rel.pred.ben <- pred.model.obj$benefits$predicted.relative.benefit
  abs.pred.ben <- pred.model.obj$benefits$predicted.absolute.benefit
  
  ## calculate observed benefit and predicted benefit for each quantile group
  # initialize
  rel.obs.ben.mat           <- as.data.frame(matrix(NA_real_, ncol(quantile.groups), 5))
  colnames(rel.obs.ben.mat) <- c("quantile", "estimate", "ci.lower", "ci.upper", "stderr")
  rel.obs.ben.mat$quantile  <- gsub(" .*$", "", colnames(quantile.groups))
  abs.obs.ben.mat  <- rel.obs.ben.mat
  abs.pred.ben.mat <- rel.obs.ben.mat
  rel.pred.ben.mat <- rel.obs.ben.mat
  or.mat           <- rel.obs.ben.mat
  
  for(i in 1:ncol(quantile.groups)){
    group <- quantile.groups[,i]
    
    # absolute observed benefit
    abs.obs.ben.mat[i, 2:5] <- 
      absolute.observed.benefit(y[group], w[group], significance.level = significance.level)
    
    # absolute predicted benefit
    abs.pred.ben.mat[i, 2:5] <- 
      predicted.benefit.inference(abs.pred.ben[group], significance.level = significance.level)
    
    # relative observed benefit
    rel.obs.ben.mat[i, 2:5] <- 
      relative.observed.benefit(y[group], w[group], significance.level = significance.level)
    
    # relative predicted benefit
    rel.pred.ben.mat[i, 2:5] <- 
      predicted.benefit.inference(rel.pred.ben[group], significance.level = significance.level)
    
    # odds ratio
    or.mat[i, 2:5] <- 
      odds.ratio(y[group], w[group], significance.level = significance.level)
    
  } # FOR
  
  return(list(absolute.observed.benefit = abs.obs.ben.mat, 
              absolute.predicted.benefit = abs.pred.ben.mat, 
              relative.observed.benefit = rel.obs.ben.mat, 
              relative.predicted.benefit = rel.pred.ben.mat, 
              odds.ratio = or.mat,
              significance.level = significance.level,
              quantile.cutoff.points = cutoffs,
              group.membership = quantile.groups))
}


#' Calculates group-level benefits from a GRF model, and the associated confidence intervals.
#' The returned benefits are the observed and predicted relative and absolute benefits as well as the odds ratio
#' 
#' @param grf.model.obj GRF model object
#' @param cutoffs the quantile cutoff points. Default is c(0.25, 0.5, 0.75), which yields the quartiles.
#' @param significance.level the significance level. Default is 0.05.
#' 
#' @export
get.benefits.grf <- function(grf.model.obj, 
                             cutoffs = c(0.25, 0.5, 0.75),
                             significance.level = 0.05){
  # calculates the observed and predicted absolute benefits along the associated confidence intervals for the GRF (relative effects cannot be computed)
  
  # extract outcome and treatment status
  y <- grf.model.obj$inputs$y.prediction.timeframe
  w <- grf.model.obj$inputs$w
  
  # group observations by their quantile of predicted baseline risk
  quantile.groups <- quantile.group(grf.model.obj$risk$risk.baseline, cutoffs)
  
  ## calculate observed benefit and predicted benefit for each quantile group
  # initialize
  abs.obs.ben.mat           <- as.data.frame(matrix(NA_real_, ncol(quantile.groups), 5))
  colnames(abs.obs.ben.mat) <- c("quantile", "estimate", "ci.lower", "ci.upper", "stderr")
  abs.obs.ben.mat$quantile  <- gsub(" .*$", "", colnames(quantile.groups))
  abs.pred.ben.mat <- abs.obs.ben.mat
  
  for(i in 1:ncol(quantile.groups)){
    group <- quantile.groups[,i]
    
    ## absolute observed benefit
    abs.obs.ben.mat[i, 2:5] <- 
      absolute.observed.benefit(y[group], w[group], significance.level = significance.level)
    
    ## absolute predicted benefit
    ate.group <- grf::average_treatment_effect(grf.model.obj$causal.forest.obj, 
                                               subset = group)
    ate <- unname(ate.group["estimate"])
    se  <- unname(ate.group["std.err"]) 
    
    # quantile of the standard normal distribution
    z <- qnorm(1 - significance.level/2)
    
    abs.pred.ben.mat[i, "estimate"] <- ate
    abs.pred.ben.mat[i, "stderr"]   <- se
    abs.pred.ben.mat[i, "ci.lower"] <- ate - z * se
    abs.pred.ben.mat[i, "ci.upper"] <- ate + z * se
    
  } # FOR
  
  return(list(absolute.observed.benefit = abs.obs.ben.mat, 
              absolute.predicted.benefit = abs.pred.ben.mat, 
              significance.level = significance.level,
              quantile.cutoff.points = cutoffs,
              group.membership = quantile.groups))
}