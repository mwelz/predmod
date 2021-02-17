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
#' ! Here we only evaluate the predictive models (effect and risk modeling).
#' ! All other models are simulated in other scripts.
#' 
#' Author: mwelz, kth
#' Last changed: February 17, 2021
#' --------------------------------------------
rm(list = ls()) ; gc() ; cat("\014")

### 0. setup ----
# load the helper functions and seeds
source(paste0(getwd(), "/funs/estimation-funs.R"))
source(paste0(getwd(), "/simulations/baseline/dgp/baseline-sim-funs.R")) 
source(paste0(getwd(), "/funs/Poiss_test_estimation-funs.R")) 
source(paste0(getwd(), "//funs/Poiss_treat_test_estimation-funs.R")) 
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

for(r in 1:R){
  
  # generate data
  data <- dgp.baseline(n = n, p = p, theta = theta, 
                       risk.rel = risk.rel, base.lifeyears=base.lifeyears, lifeyear.coefs=lifeyear.coefs, lc.effect.ly=lc.effect.ly, seed = seeds[r])
  
  ## 1.1. risk modeling ----
  # risk modeling with Lasso
  risk.model <- risk.modeling(X = data$x, w = data$w, y = data$y, alpha = 1, offset.lp = TRUE)
  
  # group-level evaluation
  rm.calib.rel <- calibration.plot(risk.model, relative = TRUE)$data
  rm.calib.abs <- calibration.plot(risk.model, relative = FALSE)$data
  
  ## 1.1.b risk modeling (poisson)----
  # risk modeling with Lasso
  risk.model.poiss <- risk.modeling.poiss(X = data$x, w = data$w, y = data$y, ly=data$ly, alpha = 1, offset.lp = TRUE)
  
  # group-level evaluation
  rm.poiss.calib.rel <- calibration.plot(risk.model.poiss, relative = TRUE)$data
  rm.poiss.calib.abs <- calibration.plot(risk.model.poiss, relative = FALSE)$data
  
  # generate arrays to store data
  if(r == 1){
    
    # sample-level
    rm.mat            <- matrix(NA_real_, R, 5)
    colnames(rm.mat)  <- c("abs.pb", "abs.pb.err", "rel.pb", "rel.pb.err", "c.index")
    em.mat            <- rm.mat
    rm.poiss.mat       <- rm.mat
    em.poiss.mat       <- rm.mat
    
    # group-level
    rm.arr            <- array(NA_real_, dim = c(4, nrow(rm.calib.rel), R))
    rownames(rm.arr)  <- c("abs.pb", "abs.ob", "rel.pb", "rel.ob")
    colnames(rm.arr)  <- as.character(rm.calib.rel$risk.quantile)
    em.arr            <- rm.arr
    rm.poiss.arr       <- rm.arr
    em.poiss.arr       <- rm.arr
    
    # variable selection in effect modeling
    em.vars           <- matrix(NA, R, 2*p+1)
    x.nam.full        <- c(paste0("x", 1:p), "w", paste0("w.x", 1:p))
    colnames(em.vars) <- x.nam.full
    
    em.poiss.vars           <- matrix(NA, R, 2*p+1)
    x.nam.full        <- c(paste0("x", 1:p), "w", paste0("w.x", 1:p))
    colnames(em.poiss.vars) <- x.nam.full
    
    # rate-ratio results
    rr.arr            <- array(NA_real_, dim = c(R,5,subgroups+1))
    colnames(rr.arr )  <- c("rateratio", "l.ci", "u.ci", "ratedifference", "p.value")
    

    
    
  } # IF
  
  # sample-level estimated relative and absolute benefit (all in absolute values)
  rm.mat[r, "abs.pb"]     <- risk.model$ate.hat
  rm.mat[r, "abs.pb.err"] <- abs(abs(data$ate) - risk.model$ate.hat)
  rm.mat[r, "rel.pb"]     <- mean(risk.model$predicted.relative.benefit)
  rm.mat[r, "rel.pb.err"] <- abs(risk.rel - mean(risk.model$predicted.relative.benefit))
  rm.mat[r, "c.index"]    <- risk.model$c.index
  
  # group-level calibration performance (in absolute values)
  rm.arr["abs.pb",,r]     <- rm.calib.abs$pb.means
  rm.arr["abs.ob",,r]     <- rm.calib.abs$ob.means
  rm.arr["rel.pb",,r]     <- rm.calib.rel$pb.means
  rm.arr["rel.ob",,r]     <- rm.calib.rel$ob.means
  
  
  # sample-level estimated relative and absolute benefit (all in absolute values)
  rm.poiss.mat[r, "abs.pb"]     <- risk.model.poiss$ate.hat
  rm.poiss.mat[r, "abs.pb.err"] <- abs(abs(data$ate) - risk.model.poiss$ate.hat)
  rm.poiss.mat[r, "rel.pb"]     <- mean(risk.model.poiss$predicted.relative.benefit)
  rm.poiss.mat[r, "rel.pb.err"] <- abs(risk.rel - mean(risk.model.poiss$predicted.relative.benefit))
  rm.poiss.mat[r, "c.index"]    <- risk.model.poiss$c.index
  
  # group-level calibration performance (in absolute values)
  rm.poiss.arr["abs.pb",,r]     <- rm.poiss.calib.abs$pb.means
  rm.poiss.arr["abs.ob",,r]     <- rm.poiss.calib.abs$ob.means
  rm.poiss.arr["rel.pb",,r]     <- rm.poiss.calib.rel$pb.means
  rm.poiss.arr["rel.ob",,r]     <- rm.poiss.calib.rel$ob.means
  
  
  
  
  
  
  ## 1.2. effect modeling ----
  # effect modeling with Lasso
  effect.model <- effect.modeling(x = data$x, w = data$w, y = data$y, alpha = 1)
  
  # sample-level estimated relative and absolute benefit (all in absolute values)
  em.mat[r, "abs.pb"]     <- effect.model$ate.hat
  em.mat[r, "abs.pb.err"] <- abs(abs(data$ate) - effect.model$ate.hat)
  em.mat[r, "rel.pb"]     <- mean(effect.model$predicted.relative.benefit)
  em.mat[r, "rel.pb.err"] <- abs(risk.rel - mean(effect.model$predicted.relative.benefit))
  em.mat[r, "c.index"]    <- effect.model$c.index
  
  # group-level calibration performance (in absolute values)
  em.calib.rel            <- calibration.plot(effect.model, relative = TRUE)$data
  em.calib.abs            <- calibration.plot(effect.model, relative = FALSE)$data
  em.arr["abs.pb",,r]     <- em.calib.abs$pb.means
  em.arr["abs.ob",,r]     <- em.calib.abs$ob.means
  em.arr["rel.pb",,r]     <- em.calib.rel$pb.means
  em.arr["rel.ob",,r]     <- em.calib.rel$ob.means
  
  # variable selection performance
  em.vars[r,]             <- x.nam.full %in% rownames(effect.model$effect.model$summary)
  
  ## 1.2.b poisson effect modeling ----
  # effect modeling with Lasso
  effect.model.poiss <- effect.modeling.poiss(X = data$x, w = data$w, y = data$y,ly=data$ly,  alpha = 1)
  
  # sample-level estimated relative and absolute benefit (all in absolute values)
  em.poiss.mat[r, "abs.pb"]     <- effect.model.poiss$ate.hat
  em.poiss.mat[r, "abs.pb.err"] <- abs(abs(data$ate) - effect.model.poiss$ate.hat)
  em.poiss.mat[r, "rel.pb"]     <- mean(effect.model.poiss$predicted.relative.benefit)
  em.poiss.mat[r, "rel.pb.err"] <- abs(risk.rel - mean(effect.model.poiss$predicted.relative.benefit))
  em.poiss.mat[r, "c.index"]    <- effect.model.poiss$c.index
  
  # group-level calibration performance (in absolute values)
  em.poiss.calib.rel            <- calibration.plot(effect.model.poiss, relative = TRUE)$data
  em.poiss.calib.abs            <- calibration.plot(effect.model.poiss, relative = FALSE)$data
  em.poiss.arr["abs.pb",,r]     <- em.calib.abs$pb.means
  em.poiss.arr["abs.ob",,r]     <- em.calib.abs$ob.means
  em.poiss.arr["rel.pb",,r]     <- em.calib.rel$pb.means
  em.poiss.arr["rel.ob",,r]     <- em.calib.rel$ob.means
  
  # variable selection performance
  em.poiss.vars[r,]             <- x.nam.full %in% rownames(effect.model.poiss$effect.model.poiss$summary)
  
 
  ## 1.3. rate-ratios ----
  
  #General rate-ratios
  rateratio.estimate              <- rate.ratio(y=data$y,w=data$w,lifeyears=data$ly)
  rr.arr[r,"rateratio",1]         <- rateratio.estimate$rate.ratio
  rr.arr[r,"l.ci",1]              <- rateratio.estimate$rate.ratio.obj$conf.int[1]
  rr.arr[r,"u.ci",1]              <- rateratio.estimate$rate.ratio.obj$conf.int[2]
  rr.arr[r,"ratedifference",1]    <- rateratio.estimate$rate.ratio.obj$estimate[2]-rateratio.estimate$rate.ratio.obj$estimate[3] #difference in rates
  rr.arr[r,"p.value",1]           <- rateratio.estimate$rate.ratio.obj$p.value
  
  #placeholder to consider immediate consideration of subgroups
#  for(a in 1:subgroups)
#  {
#    rate.ratio(y=data$y,w=data$w,lifeyears=data$ly,subgroup=data$subgroup[subgroups])
#    rr.arr[r,"rateratio",1+a]         <- rateratio.estimate$rate.ratio
#    rr.arr[r,"l.ci",1+a]              <- rateratio.estimate$rate.ratio.obj$conf.int[1]
#    rr.arr[r,"u.ci",1+a]              <- rateratio.estimate$rate.ratio.obj$conf.int[2]
#    rr.arr[r,"ratedifference",1+a]    <- rateratio.estimate$rate.ratio.obj$estimate[2]-rateratio.estimate$rate.ratio.obj$estimate[3] #difference in rates
#    rr.arr[r,"p.value",1+a]           <- rateratio.estimate$rate.ratio.obj$p.value
#    
# }
 
} # FOR

save(list = c("rm.mat", "rm.arr", "em.mat", "em.arr", "em.vars", "R", "n", "p", "theta", "risk.rel","rr.arr","rm.poiss.mat", "rm.poiss.arr", "em.poiss.mat", "em.poiss.arr", "em.poiss.vars"), 
     file = paste0(getwd(), "/simulations/baseline/sim-predictive-models/sim-results-raw.Rdata"))

# TODO: relative benefit as 1 - relative benefit to make it a risk reduction?
# TODO: 1-variable-at-a-time; survival models, and GRF. They will be done in different scripts
