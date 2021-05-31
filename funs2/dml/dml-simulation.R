rm(list = ls()) ; cat("\014")

logistic <- function(x) 1 / ( 1 + exp(-x))

# load dml function
source(paste0(getwd(), "/funs/dml/dml-funs.R"))

# Simulate data
set.seed(124341)
n <- 5000
p <- 5
R <- 100

coverage <- rep(NA, R)

for(r in 1:R){
  
  # treatment assignment (assume RCT)
  w <- rbinom(n, 1, 0.5) 
  
  # covariates
  X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
  
  # coefficients (including an intercept)
  theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)
  
  # compute Pr(Y = 1 | X) for each individual (with noise)
  pi0 <- logistic(as.numeric(cbind(1, X) %*% theta))
  
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
  
  # estimate model 
  dml.obj <-  dml(X = X, w = w, y = y, ml_m = "glm", ml_g = "glm", significance.level = 0.05)
  coverage[r] <- dml.obj$confidence.interval[1] <= ate & ate <= dml.obj$confidence.interval[2]
  
}

# average coverage (95% expected)
mean(coverage) # 0.96
