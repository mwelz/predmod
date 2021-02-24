#' --------------------------------------------
#' sim-grf.R
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
#' Last changed: February 24, 2021
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
base.lifeyears = 10 
lifeyear.coefs = c(-0.5, 0.3, -0.7, 0.1, -0.4) 
lc.effect.ly = 0.2
subgroups = 0 #placeholder value for subgroup evaluation for rate-ratios
r=1

# list with quantile groups to consider for the calibration analysis
quantile.groups.ls <- list(groups.2  = 0.5,
                           groups.4  = c(0.25, 0.5, 0.75),
                           groups.6  = c(0.1, 0.3, 0.5, 0.7, 0.9),
                           groups.8  = c(0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9),
                           groups.10 = seq(0.1, 0.9, 0.1))
quantile.groups <- names(quantile.groups.ls)

for(r in 1:R){
  
  # generate data
  data <- dgp.baseline(n = n, p = p, theta = theta, 
                       risk.rel = risk.rel, 
                       base.lifeyears = base.lifeyears, 
                       lifeyear.coefs = lifeyear.coefs,
                       lc.effect.ly = lc.effect.ly, 
                       seed = seeds[r])
  
  ## 1. GRF modeling ----
  grf.model <- grf.modeling(X = data$x, w = data$w, y = data$y, num.trees = 2000)
  
  
  # generate arrays to store data
  if(r == 1){
    
    # sample-level
    grf.mat            <- matrix(NA_real_, R, 4)
    colnames(grf.mat)  <- c("abs.pb", "abs.pb.err", "in95%ci", "c.index")

    # initialize arrays that store the group-level results 
    grf.arr.ls <- lapply(quantile.groups, function(group){
      
      grf.arr.group <<- array(NA_real_, dim = c(3, length(quantile.groups.ls[[group]]) + 1, R))
      rownames(grf.arr.group) <<- c("abs.pb", "abs.ob", "abs.cover") 
      colnames(grf.arr.group) <<- 
        gsub(" .*$", "", colnames(quantile.group(1:100, 
                                                 cutoffs = quantile.groups.ls[[group]], 
                                                 quantile.nam = TRUE)))
      grf.arr.group
    })
    names(grf.arr.ls) <- quantile.groups
    
  } # IF
  
  ### 2. evaluate sample-level results ----
  # sample-level estimated relative and absolute benefit 
  ate.ci.lo                <- grf.model$ate.hat - grf.model$ate.hat.se * qt(0.975, df = n - p)
  ate.ci.up                <- grf.model$ate.hat + grf.model$ate.hat.se * qt(0.975, df = n - p)
  grf.mat[r, "abs.pb"]     <- abs(grf.model$ate.hat)
  grf.mat[r, "abs.pb.err"] <- abs(abs(data$ate) - abs(grf.model$ate.hat))
  grf.mat[r, "in95%ci"]    <- 1 * ((ate.ci.lo <= data$ate) & (data$ate <= ate.ci.up)) # is ATE in 95% CI?
  grf.mat[r, "c.index"]    <- grf.model$c.index
  
  
  ### 3. evaluate group-level results ----
  for(group in quantile.groups){
    
    ## 3.1 risk modeling, group-level ----
    grf.calib <- calibration.plot.grf(grf.model, 
                                      quantiles = quantile.groups.ls[[group]], 
                                      alpha.significance = 0.05)$data
    
    # group-level calibration performance (in absolute values) 
    grf.arr.full          <- grf.arr.ls[[group]]
    grf.arr               <- grf.arr.full[,,r]
    grf.arr["abs.pb",]    <- grf.calib$pb.means
    grf.arr["abs.ob",]    <- grf.calib$ob.means

    # is the predicted benefit contained in CI of absolute benefit? TODO: something is weird here, especially for relative! This has to do with fact that in absolute, the CI's are in absolute value, but this is not true for Cis of relative!
    grf.arr["abs.cover",] <- 
      (grf.calib$ob.means.ci.lo <= grf.calib$pb.means) & 
      (grf.calib$pb.means <= grf.calib$ob.means.ci.up)
    
    # re-assign to storage list
    grf.arr.full[,,r]    <- grf.arr
    grf.arr.ls[[group]]  <- grf.arr.full
    
  } # FOR quantile.groups

} # FOR

save(list = c("grf.mat", "grf.arr.ls", "R", "n", "p", "theta", "risk.rel", "quantile.groups.ls"), 
     file = paste0(getwd(), "/simulations/baseline/sim-grf/sim-results-raw.Rdata"))
