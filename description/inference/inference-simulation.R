#' --------------------------------------------------------------------------------------------
#' This is a simulation study on the accuracy of confidence intervals for "observed benefit"-type estimators. The confidence intervals are found to indeed be accurate.
#'
#' Literature:
#' 
#' On the relation between risk ratio and odds ratio:
#' https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Mantel-Haenszel_Test.pdf
#' 
#' Statistical properties of the risk ratio and the odds ratio:
#' http://www.stat.ucla.edu/~dinov/courses_students.dir/05/Fall/STAT13.2.dir/STAT13_notes.dir/lecture10.pdf
#' 
#' Highly-cited article on odds ratio:
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2938757/
#' 
#' On confidence intervals for the relative risk:
#' https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_confidence_intervals/bs704_confidence_intervals8.html
#' 
#' Author: mwelz
#' Last changed: Mar 30, 2021
#' --------------------------------------------------------------------------------------------
rm(list = ls()) ; cat("\014")

# define the theoretical probabilities (the "truth")
p1 <- 0.4 # Pr(Y=1 | W=1)
p0 <- 0.7 # Pr(Y=1 | W=0)

# define the theoretical odds ratio, risk ratio, and average treatment effect
or  <- (p1 / (1 - p1) ) / (p0 / (1 - p0)) # odds ratio
rr  <- p1 / p0                            # relative risk
ate <- p1 - p0                            # average treatment effect

# initialize
set.seed(14541)
n <- 1000 # number of sampled individuals
R <- 1000 # number of simulation runs
z <- qnorm(0.975) # 95% standard normal quantile (we thus expect 95% coverage)
coverage.or <- rep(NA, R)
coverage.rr <- rep(NA, R)
coverage.ate <- rep(NA, R)

for(r in 1:R){
  
  # sample treatment status and outcome from Bernoulli distribution
  w <- rbinom(n, 1, 0.5)
  y <- rep(NA, n)
  y[w == 1] <- rbinom(sum(w), 1, p1)
  y[w == 0] <- rbinom(n-sum(w), 1, p0)
  
  # contingency table
  a <- sum(w == 1 & y == 1)
  b <- sum(w == 1 & y == 0)
  c <- sum(w == 0 & y == 1)
  d <- sum(w == 0 & y == 0)
  
  tab <- matrix(c(a,b,c,d), 2, byrow = TRUE)
  rownames(tab) <- c("w=1", "w=0")
  colnames(tab) <- c("y=1", "y=0")
  
  # estimates for the probabilities
  p1.hat <- a / (a + b)
  p0.hat <- c / (c + d)
  
  
  ## inference on the odds ratio ----
  # empirical odds ratio and its standard error
  or.hat    <- (p1.hat / (1 - p1.hat) ) / (p0.hat / (1 - p0.hat))
  se.or.hat <- sqrt(1/a + 1/b + 1/c + 1/d)
  
  # construct a confidence interval for the odds ratio
  ci.lo.or.hat <- or.hat * exp(-z * se.or.hat)
  ci.up.or.hat <- or.hat * exp(z * se.or.hat)
  
  # assess the coverage
  coverage.or[r] <- (ci.lo.or.hat <= or) & (or <= ci.up.or.hat)
  
  
  ## inference on the risk ratio ----
  # empirical risk ratio and its standard error
  rr.hat    <- p1.hat / p0.hat
  se.rr.hat <- sqrt(b / (a * (a + b)) + d / (c * (c + d)))
  
  # construct a confidence interval for the risk ratio
  ci.lo.rr.hat <- rr.hat * exp(-z * se.rr.hat)
  ci.up.rr.hat <- rr.hat * exp(z * se.rr.hat)
  
  # assess the coverage
  coverage.rr[r] <- (ci.lo.rr.hat <= rr) & (rr <= ci.up.rr.hat)
  
  
  ## average treatment effect ----
  # use the t.test() function to calculate Welch standard error
  ttest      <- t.test(y[w==1], y[w==0], var.equal = FALSE) # TODO: add 'false' to main code. 
  ate.hat    <- unname(ttest$estimate[1] - ttest$estimate[2])
  se.ate.hat <- unname(ttest$stderr) 
  
  # 95% quantile of t-distribution
  t <- qt(0.975, df = ttest$parameter)
  
  # construct confidence interval for the average treatment effect
  ci.lo.ate.hat <- ate.hat - t * se.ate.hat
  ci.up.ate.hat <- ate.hat + t * se.ate.hat
  
  # assess the coverage
  coverage.ate[r] <- (ci.lo.ate.hat <= ate) & (ate <= ci.up.ate.hat)

}

# evaluation: we expect approx. 95% coverage everywhere below
mean(coverage.or)  # 0.948, so the CI for odds ratio is indeed accurate
mean(coverage.rr)  # 0.956, so the CI for risk ratio is indeed accurate
mean(coverage.ate) # 0.953, so the CI for ATE is indeed accurate