#' baseline DGP
dgp.baseline <- function(n, p, theta = -2.5, 
                         time.at.risk.baseline = 7,
                         survivors.additional.time = 3){
  
  # sample covariates
  X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
  
  # sample treatment asisgnment
  W <- rbinom(n, size = 1, prob = 0.5)
  
  # Rademacher-type coefficients of covariates
  beta <- rep(-1, p)
  beta[1:p %% 2 == 0] <- 1
  beta0 <- 0
  
  # predictive function
  eta0 <- beta0 + as.numeric(X %*% beta)
  eta1 <- eta0 + theta
  
  # mortality risk
  p0 <- plogis(eta0)
  p1 <- plogis(eta1)
  
  # potential outcomes
  Y0 <- rbinom(n, size = 1, prob = p0)
  Y1 <- rbinom(n, size = 1, prob = p1)
  
  # observed outcome
  Y <- ifelse(W == 1, Y1, Y0)
  
  # time at risk
  T. <- rep(time.at.risk.baseline, n) + ifelse(Y == 0, survivors.additional.time, 0)
 
  # return
  return(list(X = X, Y = Y, W = W, time = T., 
              ATE = mean(p1 - p0),
              ARTE = mean(p1 / p0),
              p1 = p1, p0 = p0,
              Y1 = Y1, Y0 = Y0))
  
} # FUN


seedmaker <- function(R) round(1e7 * runif(R), 0)