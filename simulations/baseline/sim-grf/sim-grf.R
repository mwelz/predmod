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
#' ! Here we only evaluate the GRF (can only model absolute benefit).
#' ! All other models are simulated in other scripts.
#' 
#' Author: mwelz
#' Last changed: February 11, 2021
#' --------------------------------------------
rm(list = ls()) ; gc() ; cat("\014")

### 0. setup ----
# load the helper functions and seeds
source(paste0(getwd(), "/funs/estimation-funs.R"))
source(paste0(getwd(), "/simulations/baseline/dgp/baseline-sim-funs.R")) 
load(paste0(getwd(), "/simulations/baseline/dgp/seeds.Rdata")) # load seeds

### 1. simulation ----
# parameters
R        <- 100
n        <- 10000
p        <- 5
theta    <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)
risk.rel <- 0.7 # relative risk (i.e. pi0 / pi1), so not the reduction!

for(r in 1:R){
  
  # generate data
  data <- dgp.baseline(n = n, p = p, theta = theta, 
                       risk.rel = risk.rel, seed = seeds[r])
  
  ## 1. GRF modeling ----
  grf.model <- grf.modeling(X = data$x, w = data$w, y = data$y, num.trees = 2000)
  
  # group-level evaluation
  grf.calib <- calibration.plot.grf(grf.model)$data
  
  # generate arrays to store data
  if(r == 1){
    
    # sample-level
    grf.mat            <- matrix(NA_real_, R, 4)
    colnames(grf.mat)  <- c("abs.pb", "abs.pb.err", "in95%ci", "c.index")

    # group-level
    grf.arr            <- array(NA_real_, dim = c(2, nrow(grf.calib), R))
    rownames(grf.arr)  <- c("abs.pb", "abs.ob")
    colnames(grf.arr)  <- as.character(grf.calib$risk.quantile)
  } # IF
  
  # sample-level estimated relative and absolute benefit (all in absolute values)
  ate.ci.lo                <- grf.model$ate.hat - grf.model$ate.hat.se * qt(0.975, df = n - p)
  ate.ci.up                <- grf.model$ate.hat + grf.model$ate.hat.se * qt(0.975, df = n - p)
  grf.mat[r, "abs.pb"]     <- abs(grf.model$ate.hat)
  grf.mat[r, "abs.pb.err"] <- abs(abs(data$ate) - abs(grf.model$ate.hat))
  grf.mat[r, "in95%ci"]    <- 1 * ((ate.ci.lo <= data$ate) & (data$ate <= ate.ci.up)) # is ATE in 95% CI?
  grf.mat[r, "c.index"]    <- grf.model$c.index
  
  # group-level calibration performance (in absolute values)
  grf.arr["abs.pb",,r]     <- grf.calib$pb.means
  grf.arr["abs.ob",,r]     <- grf.calib$ob.means

} # FOR

save(list = c("grf.mat", "grf.arr", "R", "n", "p", "theta", "risk.rel"), 
     file = paste0(getwd(), "/simulations/baseline/sim-grf/sim-results-raw.Rdata"))

# TODO: relative benefit as 1 - relative benefit to make it a risk reduction?
# TODO: 1-variable-at-a-time; survival models, and GRF. They will be done in different scripts
