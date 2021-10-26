#' Partition a vector into groups based on its quantiles.
#'
#' Partitions a vector into quantile groups and returns a logical matrix indicating group membership.
#'
#' @param x The vector to be partitioned
#' @param cutoffs The quantile cutoffs for the partition. Default are the quartiles: \code{c(0.25, 0.5, 0.75)}.
#' @param names_quantile Logical. If \code{TRUE}, then the column names of the returned matrix are the quantiles as in \code{cutoffs}. If \code{FALSE}, the names are the numeric intervals that constitute the grouping.
#'
#' @return An object of the class \code{quantile_group}, which is a logical matrix indicating group membership
#'
#' @export
quantile_group <- function(x,
                           cutoffs = c(0.25, 0.5, 0.75),
                           names_quantile = TRUE){
  # cutoffs are the quantile cutoffs (like c(0.25, 0.5, 0.75))
  
  # get quatiles
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
  group.mat <- matrix(NA, length(x), length(group.nam))
  nam       <- rep(NA, length(group.nam))
  
  for(j in 1:length(group.nam)){
    if(j == 1){
      nam[j] <- paste0("<", 100*cutoffs[j], "% quantile")
    } else if (j == length(group.nam)){
      nam[j] <- paste0(">=", 100*cutoffs[j-1], "% quantile")
    } else{
      nam[j] <- paste0("[", 100*cutoffs[j-1], ",", 100*cutoffs[j], ")% quantile")
    }
    group.mat[,j] <- groups == group.nam[j]
  }
  
  if(names_quantile){
    colnames(group.mat) <- nam
  } else{
    colnames(group.mat) <- group.nam
  }
  return(structure(group.mat, type = "quantile_group"))
} # FUN


#' calculates the absolute observed benefit, which corresponds to the difference in \code{mean(y[W=w])}, and the associated confidence interval.
#' 
#' @param y vector of binary outcomes
#' @param w vector of binary treatment status
#' @param significance.level the significance level. Default is 0.05.
#'
#' @export
absolute.observed.benefit <- function(y, w, significance.level = 0.05){
  
  ttest  <- stats::t.test(y[w==1], y[w==0], paired = FALSE, var.equal = FALSE) 
  aob    <- unname(ttest$estimate[1] - ttest$estimate[2])
  se.aob <- unname(ttest$stderr) 
  
  # quantile of the standard normal distribution
  z <- stats::qnorm(1 - significance.level/2)
  
  return(c(estimate = aob, 
           ci.lower = aob - z * se.aob,
           ci.upper = aob + z * se.aob,
           stderr   = se.aob))
} # FUN


#' calculates the relative observed benefit, which corresponds to a risk ratio, and the associated confidence interval
#' 
#' @param y vector of binary outcomes
#' @param w vector of binary treatment status
#' @param significance.level the significance level. Default is 0.05.
#'
#' @export
relative.observed.benefit <- function(y, w, significance.level = 0.05){
  
  # contingency table with the following structure
  #      Y=1 Y=0
  #  W=1  a   b
  #  W=0  c   d
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
  z <- stats::qnorm(1 - significance.level/2)
  
  return(c(estimate = rob, 
           ci.lower = rob * exp(-z * se.rob),
           ci.upper = rob * exp(z * se.rob),
           stderr   = se.rob))
} # FUN



#' calculates the odds ratio, and the associated confidence interval. 
#' 
#' @param y vector of binary outcomes
#' @param w vector of binary treatment status
#' @param significance.level the significance level. Default is 0.05.
#'
#' @export
odds.ratio <- function(y, w, significance.level = 0.05){
  
  # contingency table with the following structure
  #      Y=1 Y=0
  #  W=1  a   b
  #  W=0  c   d
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
  z <- stats::qnorm(1 - significance.level/2)
  
  return(c(estimate = or, 
           ci.lower = or * exp(-z * se.or),
           ci.upper = or * exp(z * se.or),
           stderr   = se.or))
} # FUN



#' performs inference on the per-individual predicted benefits
#' 
#' @param predicted.benefits a vector of per-individual predicted benefits (either absolute of relative)
#' @param significance.level the significance level. Default is 0.05.
#' 
#' @export
predicted.benefit.inference <- function(predicted.benefits, significance.level = 0.05){
  
  ttest <- stats::t.test(predicted.benefits)
  pb    <- unname(ttest$estimate)
  se.pb <- unname(ttest$stderr)
  
  # quantile of the standard normal distribution
  z <- stats::qnorm(1 - significance.level/2)
  
  return(c(estimate = pb, 
           ci.lower = pb - z * se.pb,
           ci.upper = pb + z * se.pb,
           stderr   = se.pb))
} # FUN


#' Calculates group-level benefits from a prediction model, and the associated confidence intervals.
#' The returned benefits are the observed and predicted relative and absolute benefits as well as the odds ratio
#' 
#' @param pred.model.obj prediction model object
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
  quantile.groups <- quantile_group(pred.model.obj$risk$risk.baseline, cutoffs)
  
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
} # FUN


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
  quantile.groups <- quantile_group(grf.model.obj$risk$risk.baseline, cutoffs)
  
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
    z <- stats::qnorm(1 - significance.level/2)
    
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
#' @param significance.level the significance level. Default is 0.05.
#' 
#' @export
get.benefits_imputation.accounter <- function(pred.model.objs_imputed, 
                                              cutoffs = c(0.25, 0.5, 0.75),
                                              significance.level = 0.05){
  
  # initialize 
  m                         <- length(pred.model.objs_imputed)
  group.names               <- gsub(" .*$", "", colnames(quantile_group(seq(0, 1, by = 0.001), cutoffs = cutoffs)))
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
    
    # group observations by their quantile of predicted baseline risk
    quantile.groups <- quantile_group(pred.model.obj$risk$risk.baseline, cutoffs)
    
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
#' @param significance.level the significance level. Default is 0.05.
#' 
#' @export
get.benefits.grf_imputation.accounter <- function(grf.model.obj_imputed, 
                                                  cutoffs = c(0.25, 0.5, 0.75),
                                                  significance.level = 0.05){
  
  # initialize 
  m                         <- length(grf.model.obj_imputed)
  group.names               <- gsub(" .*$", "", colnames(quantile_group(seq(0, 1, by = 0.001), cutoffs = cutoffs)))
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
    
    # group observations by their quantile of predicted baseline risk
    quantile.groups <- quantile_group(grf.model.obj$risk$risk.baseline, cutoffs)
    
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