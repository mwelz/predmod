#' --------------------------------------------
#' baseline-dgp-sim.R
#' This is the baseline DGP on which the classical methods are expected to do well.
#' Baseline DGP: linear, no interaction effects, constant HTE in the sense that there
#' is a constant _relative_ treatment effect.
#' 
#' Author: mwelz
#' Last changed: December 15, 2020
#' --------------------------------------------

rm(list = ls()) ; gc() ; cat("\014")

### 0. setup ----
# load the helper functions
source(paste0(getwd(), "/funs/estimation-funs.R"))

# create dgp function
dgp.baseline <- function(n, p,
                         theta = c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4), # icnludes intercept
                         risk.rel = 0.7,
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


### 1. simulation ----
# parameters
R        <- 100
n        <- 10000
p        <- 5
theta    <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)
risk.rel <- 0.7 # relative risk (i.e. pi0 / pi1), so not the reduction!

# seeds
set.seed(83267)
seeds <- round(1e7 * runif(R), 0)

r=1
for(r in 1:R){
  
  # generate data
  data <- dgp.baseline(n = n, p = p, theta = theta, 
                       risk.rel = risk.rel, seed = seeds[r])
  
  ## 1.1. risk modeling ----
  # risk modeling with Lasso
  risk.model <- risk.modeling(X = data$x, w = data$w, y = data$y, alpha = 1, offset.lp = TRUE)
  
  # group-level evaluation
  rm.calib.rel <- calibration.plot(risk.model, relative = TRUE)$data
  rm.calib.abs <- calibration.plot(risk.model, relative = FALSE)$data
  
  # generate arrays to store data
  if(r == 1){
    
    # sample-level
    rm.mat            <- matrix(NA_real_, R, 4)
    colnames(rm.mat)  <- c("abs.pb", "abs.pb.err", "rel.pb", "rel.pb.err")
    em.mat            <- rm.mat
    
    # group-level
    rm.arr            <- array(NA_real_, dim = c(4, nrow(rm.calib.rel), R))
    rownames(rm.arr)  <- c("abs.pb", "abs.ob", "rel.pb", "rel.ob")
    colnames(rm.arr)  <- as.character(rm.calib.rel$risk.quantile)
    em.arr            <- rm.arr
    
    # variable selection in effect modeling
    em.vars           <- matrix(NA, R, 2*p+1)
    x.nam.full        <- c(paste0("x", 1:p), "w", paste0("w.x", 1:p))
    colnames(em.vars) <- x.nam.full
    
  } # IF
  
  # sample-level estimated relative and absolute benefit (all in absolute values)
  rm.mat[r, "abs.pb"]     <- risk.model$ate.hat
  rm.mat[r, "abs.pb.err"] <- abs(abs(data$ate) - risk.model$ate.hat)
  rm.mat[r, "rel.pb"]     <- mean(risk.model$predicted.relative.benefit)
  rm.mat[r, "rel.pb.err"] <- abs(risk.rel - mean(risk.model$predicted.relative.benefit))
  
  # group-level calibration performance (in absolute values)
  rm.arr["abs.pb",,r]     <- rm.calib.abs$pb.means
  rm.arr["abs.ob",,r]     <- rm.calib.abs$ob.means
  rm.arr["rel.pb",,r]     <- rm.calib.rel$pb.means
  rm.arr["rel.ob",,r]     <- rm.calib.rel$ob.means
  
  
  ## 1.2. effect modeling ----
  # effect modeling with Lasso
  effect.model <- effect.modeling(x = data$x, w = data$w, y = data$y, alpha = 1)
  
  # sample-level estimated relative and absolute benefit (all in absolute values)
  em.mat[r, "abs.pb"]     <- effect.model$ate.hat
  em.mat[r, "abs.pb.err"] <- abs(abs(data$ate) - effect.model$ate.hat)
  em.mat[r, "rel.pb"]     <- mean(effect.model$predicted.relative.benefit)
  em.mat[r, "rel.pb.err"] <- abs(risk.rel - mean(effect.model$predicted.relative.benefit))
  
  # group-level calibration performance (in absolute values)
  em.calib.rel            <- calibration.plot(effect.model, relative = TRUE)$data
  em.calib.abs            <- calibration.plot(effect.model, relative = FALSE)$data
  em.arr["abs.pb",,r]     <- em.calib.rel$pb.means
  em.arr["abs.ob",,r]     <- em.calib.rel$ob.means
  em.arr["rel.pb",,r]     <- em.calib.rel$pb.means
  em.arr["rel.ob",,r]     <- em.calib.rel$ob.means
  
  # variable selection performance
  em.vars[r,]             <- x.nam.full %in% rownames(effect.model$effect.model$summary)
  
  # TODO: relative benefit as 1 - relative benefit to make it a risk reduction?
  # TODO: 1-variable-at-a-time; survival models, and GRF
  
} # FOR

