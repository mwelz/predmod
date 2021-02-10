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
#' Author: mwelz
#' Last changed: February 10, 2021
#' --------------------------------------------
dgp.baseline <- function(n, p,
                         theta = c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4), # includes intercept
                         risk.rel = 0.7,  # relative risk (i.e. pi0 / pi1), so not the reduction!
                         seed = 1){
  
  set.seed(seed)
  
  # treatment assignment (assume RCT)
  w <- rbinom(n, 1, 0.5) 
  
  # covariates
  x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
  
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
  
  return(list(x = x, w = w, y = y, y0 = y0, y1 = y1, 
              ate = ate, hte = hte, pi0 = pi0, pi1 = pi1))
}
