rm(list = ls()) ; cat("\014")

get.benefits.grf <- function(grf.model.obj, 
                         cutoffs = c(0.25, 0.5, 0.75)){
  # cannot give relative risk!
  y <- grf.model.obj$inputs$y
  w <- grf.model.obj$inputs$w
  
  # group observations by their quantile of predicted baseline risk (predictions from random forest)
  quantile.groups <- quantile.group(grf.model.obj$risk.baseline, cutoffs)
  
  # get predicted absolute benefit (predictions from causal forest)
  pred.ben <- grf.model.obj$predicted.absolute.benefit
  
  ## calculate observed benefit and predicted benefit for each quantile group
  # initialize
  obs.ben.mat           <- as.data.frame(matrix(NA_real_, ncol(quantile.groups), 4))
  colnames(obs.ben.mat) <- c("quantile", "mean", "stderr", "df")
  obs.ben.mat$quantile  <- gsub(" .*$", "", colnames(quantile.groups))
  pred.ben.mat          <- obs.ben.mat
  
  for(i in 1:ncol(quantile.groups)){
    group <- quantile.groups[,i]
    
    ## observed absolute benefit (grf cannot do relative benefit)
    # corresponds to difference in mean(y[group & W=w])
    ttest                    <- t.test(y[group & w == 1], y[group & w == 0])
    obs.ben.mat[i, "mean"]   <- unname(ttest$estimate[1] - ttest$estimate[2])
    obs.ben.mat[i, "stderr"] <- unname(ttest$stderr) 
    obs.ben.mat[i, "df"]     <- unname(ttest$parameter)
    
    ## group by predicted benefit
    ate.group <- grf::average_treatment_effect(grf.model.obj$causal.forest.obj, 
                                               subset = group)
    pred.ben.mat[i, "mean"]   <- unname(ate.group["estimate"])
    pred.ben.mat[i, "stderr"] <- unname(ate.group["std.err"]) 
    pred.ben.mat[i, "df"]     <- sum(group) - ncol(grf.model.obj$inputs$X)
  } # IF
  
  return(list(
    group.observed.benefit = obs.ben.mat,
    group.predicted.benefit = pred.ben.mat,
    group.membership = quantile.groups
  ))
}


calibration.plot.grf <- function(grf.model.obj,
                             quantiles = c(0.25, 0.5, 0.75), 
                             title = NULL,
                             alpha.significance = 0.05){
  # only absolute benefit possible, not relative!
  # get observed and predicted benefit by quantile group
  benefits <- get.benefits.grf(grf.model.obj = grf.model.obj, 
                           cutoffs = quantiles)
  
  # make everything positive for visualization
  benefits$group.predicted.benefit$mean <- abs(benefits$group.predicted.benefit$mean) 
  benefits$group.observed.benefit$mean  <- abs(benefits$group.observed.benefit$mean)
  
  # limits of the plot
  limits <- c(min(benefits$group.predicted.benefit$mean, benefits$group.observed.benefit$mean, -0.2),
              max(benefits$group.predicted.benefit$mean, benefits$group.observed.benefit$mean, 0.3))
  
  # get the whiskers (SE*quantile)
  whisker.obs.ben <- benefits$group.observed.benefit$stderr *
    qt(1-alpha.significance/2, benefits$group.observed.benefit$df)
  
  # the plot
  library(ggplot2)
  
  # make sure risk quantile is in correct order
  risk.quantile <- factor(benefits$group.observed.benefit$quantile)
  lv <- levels(risk.quantile)
  lv <- lv[order.intervals(lv, quantile.nam = TRUE)]
  risk.quantile <- factor(risk.quantile, levels = lv)
  
  df <- data.frame(pb.means = benefits$group.predicted.benefit$mean,
                   ob.means = benefits$group.observed.benefit$mean,
                   ob.means.ci.up = benefits$group.observed.benefit$mean + whisker.obs.ben,
                   ob.means.ci.lo = benefits$group.observed.benefit$mean - whisker.obs.ben,
                   risk.quantile = risk.quantile)
  
  if(is.null(title)) title <- "Calibration plot"
  
  ggplot(mapping = aes(x = pb.means,
                       y = ob.means, color = risk.quantile), data = df) +
    geom_point() +
    geom_linerange(mapping = aes(ymin = ob.means.ci.lo,
                                 ymax = ob.means.ci.up)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = limits, ylim = limits) +
    labs(x = "Predicted absolute benefit",
         y = "Observed benefit") +
    theme_classic() +
    ggtitle(title) +
    theme(legend.position = "bottom")
}


grf.modeling <- function(X, y, w, num.trees.cf = 2000, num.trees.rf = 2000,...){
  # no relative risk modeling possible!
  # get causal forest (for predicted benefit) and random forest (for risk estimation)
  cf <- grf::causal_forest(X = X, Y = y, W = w, num.trees = num.trees.cf, ...)
  rf <- grf::regression_forest(X = X, Y = y, num.trees = num.trees.rf, ...)
  
  # initialize object
  grf.model.obj <- list()
  
  # random forest's precitions is baseline risk estimation
  grf.model.obj$risk.baseline <- as.numeric(rf$predictions)
  
  # causal forest's individual treatment effect estimates are predicted absolute benefit
  grf.model.obj$predicted.absolute.benefit <- as.numeric(cf$predictions)
  
  # collect everything relevant and return
  grf.model.obj$inputs            <- list(X = x, w = w, y = y)
  grf.model.obj$causal.forest.obj <- cf
  grf.model.obj$ate               <- unname(grf::average_treatment_effect(cf)["estimate"])
  
  # TODO: experimental: C statistic (not sure if this is correct as we are using baseline risk)
  # c.index <- unname(Hmisc::rcorr.cens(x = grf.model.obj$risk.baseline, S = y)[1])
  
  return(grf.model.obj)
}

### test ----
# load the helper functions
source(paste0(getwd(), "/funs/estimation-funs.R"))

# define the logistic function
logistic <- function(x) 1 / (1 + exp(-x))

set.seed(2)
n <- 10000
p <- 5

# treatment assignment (assume RCT)
w <- rbinom(n, 1, 0.5) 

# covariates
x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))

# coefficients (including an intercept)
theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)

# compute Pr(Y = 1 | X) for each individual (with noise)
eps <-  rnorm(n, mean = 0, sd = 0.5)
pi0 <- logistic(as.numeric(cbind(1, x) %*% theta) + eps)

# relative constant risk reduction of 30%
scaling <- 0.7
pi1 <- pi0 * scaling 
hte <- pi1 - pi0

# true ATE
ate <- mean(pi1) - mean(pi0)

# create binary outcomes
y0 <- rbinom(n, 1, pi0)
y1 <- rbinom(n, 1, pi1)
y  <- ifelse(w == 1, y1, y0) # observed outcome

## GRF test ----
grf.obj <- grf.modeling(X = x, y = y, w = w)
grf.calibration <- calibration.plot.grf(grf.obj)
