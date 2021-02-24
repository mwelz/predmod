#' --------------------------------------------
#' baseline-dgp-sim.R
#' This is the baseline DGP on which the classical methods are expected to do well.
#' Baseline DGP: linear, no interaction effects, constant HTE in the sense that there
#' is a constant _relative_ treatment effect.
#' The DGP has the following properties:
#' - it is an RCT
#' - 5 iid normally distributed covariates
#' - constant relative risk reduction of 30%. That is, 
#'   Pr(Y = 1 | X, W = 1) = 0.7 * Pr(Y = 1 | X, W = 0),
#'   where Y = 1 indicates a death.
#' 
#' 
#' 
#' Author: mwelz, kth
#' Last changed: February 17, 2021 by kth
#' --------------------------------------------
dgp.baseline <- function(n, p,
                         theta = c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4), # includes intercept
                         risk.rel = 0.7,  # relative risk (i.e. pi0 / pi1), so not the reduction!
                         base.lifeyears = 10, #Basic number of lifeyears per person before adjustment for other covariates
                         lifeyear.coefs = c(-0.5, 0.3, -0.7, 0.1, -0.4), #coefficients related to lifeyears; should match number of covariates in theta (aside from intercept). In general, a reasonable assumption would be that a higher risk for lung cancer is correlated with a lower number of remaining lifeyears (i.e. reversed sign of risk for lung cancer).
                         lc.effect.ly = 0.2, #Effect of lung cancer mortality on remaining lifeyears. We may want to replace this in the future with a survival distribution to make this more generalizable.
                         seed = 1){
  
  # NB: we could potentially add a variable  trial.period, which represents the number of years of follow-up of the trial. Individuals with lifyears >= trial.period 
  # would have their lifeyears set to trial.period, as their remaining lifeyears are censored. However, this would also require modelling a time of lung cancer death
  # which would affect the outcome variable; i.e. y = 1 if the time of lc death is <= trial.period and y = o if it is > trial.period.
  # this would be best integrated along with a survival function.
  

  set.seed(seed)
  
  # treatment assignment (assume RCT)
  w <- rbinom(n, 1, 0.5) 
  
  # covariates
  x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
  colnames(x) <- paste0("x", 1:p)
  
  # define the logistic function
  logistic <- function(x) 1 / (1 + exp(-x))
  
  # compute Pr(Y = 1 | X) for each individual (with noise)
  eps <- rnorm(n, mean = 0, sd = 0.5)
  pi0 <- logistic(as.numeric(cbind(1, x) %*% theta) + eps)
  
  # relative constant risk reduction
  pi1 <- pi0 * risk.rel 
  hte <- pi1 - pi0
  
  # true ATE
  ate <- mean(pi1) - mean(pi0)
  
  # create binary outcomes
  y0 <- rbinom(n, 1, pi0)
  y1 <- rbinom(n, 1, pi1)
  y  <- ifelse(w == 1, y1, y0) # observed outcome
  
  #Generate lifeyears
  ly_intermediate <- as.numeric((base.lifeyears+ x  %*%lifeyear.coefs) * ifelse(y == 1, lc.effect.ly, 1)) 
  # set LY of 0 or less to 0 LY
  ly <- ifelse(ly_intermediate <=0, 0, ly_intermediate) 																																												  
  
    return(list(x = x, w = w, y = y, y0 = y0, y1 = y1, 
              ate = ate, hte = hte, pi0 = pi0, pi1 = pi1, ly=ly))
}
