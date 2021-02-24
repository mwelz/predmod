#' --------------------------------------------
#' sim-predictive-models.R
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
#' ! Here we only evaluate the predictive models (effect and risk modeling).
#' ! I moved the Poisson models to a separate file (mwelz)
#' ! All other models are simulated in other scripts.
#' 
#' Author: mwelz, kth
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
  
  ### 1.1 risk & effect modeling ----
  # risk modeling with Lasso
  risk.model <- risk.modeling(X = data$x, w = data$w, y = data$y, alpha = 1, offset.lp = TRUE)
  
  # effect modeling with Lasso
  effect.model <- effect.modeling(X = data$x, w = data$w, y = data$y, alpha = 1)
  
  
  # generate arrays to store data
  if(r == 1){
    
    # sample-level
    rm.mat            <- matrix(NA_real_, R, 5)
    colnames(rm.mat)  <- c("abs.pb", "abs.pb.err", "rel.pb", "rel.pb.err", "c.index")
    em.mat            <- rm.mat

    # initialize arrays that store the group-level results 
    em.arr.ls <- lapply(quantile.groups, function(group){
      
      rm.arr.group <<- array(NA_real_, dim = c(6, length(quantile.groups.ls[[group]]) + 1, R))
      rownames(rm.arr.group) <<- c("abs.pb", "abs.ob", "abs.cover", "rel.pb", "rel.ob", "rel.cover") 
      colnames(rm.arr.group) <<- 
        gsub(" .*$", "", colnames(quantile.group(1:100, 
                                                 cutoffs = quantile.groups.ls[[group]], 
                                                 quantile.nam = TRUE)))
      rm.arr.group
    })
    names(em.arr.ls) <- quantile.groups
    rm.arr.ls        <- em.arr.ls
    
    # variable selection in effect modeling
    em.vars           <- matrix(NA, R, 2*p+1)
    x.nam.full        <- c(colnames(data$x), "w", paste0("w.", colnames(data$x)))
    colnames(em.vars) <- x.nam.full
    
  } # IF
  
  
  ### 2. evaluate sample-level results ----
  # sample-level estimated relative and absolute benefit (all in absolute values)
  rm.mat[r, "abs.pb"]     <- risk.model$ate.hat
  rm.mat[r, "abs.pb.err"] <- abs(abs(data$ate) - risk.model$ate.hat)
  rm.mat[r, "rel.pb"]     <- mean(risk.model$predicted.relative.benefit)
  rm.mat[r, "rel.pb.err"] <- abs(risk.rel - mean(risk.model$predicted.relative.benefit))
  rm.mat[r, "c.index"]    <- risk.model$c.index
  
  # sample-level estimated relative and absolute benefit (all in absolute values)
  em.mat[r, "abs.pb"]     <- effect.model$ate.hat
  em.mat[r, "abs.pb.err"] <- abs(abs(data$ate) - effect.model$ate.hat)
  em.mat[r, "rel.pb"]     <- mean(effect.model$predicted.relative.benefit)
  em.mat[r, "rel.pb.err"] <- abs(risk.rel - mean(effect.model$predicted.relative.benefit))
  em.mat[r, "c.index"]    <- effect.model$c.index
  
  # variable selection performance of effect modeling
  em.vars[r,]             <- x.nam.full %in% rownames(effect.model$effect.model$summary)
  
  
  ### 3. evaluate group-level results ----
  for(group in quantile.groups){

    ## 3.1 risk modeling, group-level ----
    rm.calib.rel <- calibration.plot(risk.model, 
                                     quantiles = quantile.groups.ls[[group]], 
                                     relative = TRUE,
                                     alpha.significance = 0.05)$data
    
    rm.calib.abs <- calibration.plot(risk.model, 
                                     quantiles = quantile.groups.ls[[group]], 
                                     relative = FALSE,
                                     alpha.significance = 0.05)$data
    
    # group-level calibration performance (in absolute values) 
    rm.arr.full          <- rm.arr.ls[[group]]
    rm.arr               <- rm.arr.full[,,r]
    rm.arr["abs.pb",]    <- abs(rm.calib.abs$pb.means)
    rm.arr["abs.ob",]    <- abs(rm.calib.abs$ob.means)
    rm.arr["rel.pb",]    <- abs(rm.calib.rel$pb.means)
    rm.arr["rel.ob",]    <- abs(rm.calib.rel$ob.means)
    
    # is the predicted benefit contained in CI of absolute benefit? TODO: something is weird here, especially for relative! This has to do with fact that in absolute, the CI's are in absolute value, but this is not true for Cis of relative!
    rm.arr["abs.cover",] <- 
    (rm.calib.abs$ob.means.ci.lo <= abs(rm.calib.abs$pb.means)) & 
      (abs(rm.calib.abs$pb.means) <= rm.calib.abs$ob.means.ci.up)
    
    rm.arr["rel.cover",] <- 
      (rm.calib.rel$ob.means.ci.lo <= rm.calib.rel$pb.means) & 
      (rm.calib.rel$pb.means <= rm.calib.rel$ob.means.ci.up)
    
    # re-assign to storage list
    rm.arr.full[,,r]    <- rm.arr
    rm.arr.ls[[group]]  <- rm.arr.full
    
    ## 3.2 effect modeling, group-level ----
    em.calib.rel <- calibration.plot(effect.model, 
                                     quantiles = quantile.groups.ls[[group]], 
                                     relative = TRUE,
                                     alpha.significance = 0.05)$data
    
    em.calib.abs <- calibration.plot(effect.model, 
                                     quantiles = quantile.groups.ls[[group]], 
                                     relative = FALSE,
                                     alpha.significance = 0.05)$data
    
    # group-level calibration performance (in absolute values) 
    em.arr.full          <- em.arr.ls[[group]]
    em.arr               <- em.arr.full[,,r]
    em.arr["abs.pb",]    <- abs(em.calib.abs$pb.means)
    em.arr["abs.ob",]    <- abs(em.calib.abs$ob.means)
    em.arr["rel.pb",]    <- abs(em.calib.rel$pb.means)
    em.arr["rel.ob",]    <- abs(em.calib.rel$ob.means)
    
    # is the predicted benefit contained in CI of absolute benefit? TODO: something is weird here, especially for relative! This has to do with fact that in absolute, the CI's are in absolute value, but this is not true for Cis of relative!
    em.arr["abs.cover",] <- 
      (em.calib.abs$ob.means.ci.lo <= abs(em.calib.abs$pb.means)) & 
      (abs(em.calib.abs$pb.means) <= em.calib.abs$ob.means.ci.up)
    
    em.arr["rel.cover",] <- 
      (em.calib.rel$ob.means.ci.lo <= em.calib.rel$pb.means) & 
      (em.calib.rel$pb.means <= em.calib.rel$ob.means.ci.up)
    
    # re-assign to storage list
    em.arr.full[,,r]    <- em.arr
    em.arr.ls[[group]]  <- em.arr.full
    
  } # FOR quantile.groups

} # FOR r


save(list = c("rm.mat", "rm.arr.ls", "em.mat", "em.arr.ls", "em.vars", "R", "n", "p", "theta", "risk.rel", "quantile.groups.ls"), 
     file = paste0(getwd(), "/simulations/baseline/sim-predictive-models/sim-results-raw.Rdata"))
