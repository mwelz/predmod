rm(list = ls()) ; cat("\014")

# load the helper functions
source(paste0(getwd(), "/risk-modeling/rm-funs/rm-funs.R"))

# define the logistic function
logistic <- function(x) 1 / (1 + exp(-x))

### 0. Data generation ----
set.seed(1)
n <- 10000
p <- 5

# treatment assignment (assume RCT)
w <- rbinom(n, 1, 0.5) 

# true treatment effect
tau <- -2

# covariates
x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))

# coefficients (including an intercept)
theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)

# compute Pr(Y = 1 | X) for each individual (with noise)
eps <-  rnorm(n, mean = 0, sd = 0.5)
pi1 <- logistic(tau + as.numeric(cbind(1, x) %*% theta) + eps)
pi0 <- logistic(as.numeric(cbind(1, x) %*% theta) + eps)
hte <- pi1 - pi0

# true ATE
ate <- mean(pi1) - mean(pi0)

# create binary outcomes
y0 <- rbinom(n, 1, pi0)
y1 <- rbinom(n, 1, pi1)
y  <- ifelse(w == 1, y1, y0) # observed outcome

# sanity check 
ate.hat <- mean(y[w==1]) - mean(y[w==0])
ate.hat - ate # accurate estimation!

### 1. risk modeling ----
risk.model <- risk.modeling(X = x, w = w, y = y, alpha = 1, offset.lp = TRUE)

# histogram of predicted benefit
hist(risk.model$predicted.benefit, main = "Histogram of Predicted Benefit", xlab = "Predicted Benefit")

# calibration plot with the 4 quartiles
calibration.plot(risk.model)

# c index (based on second stage risk probabilities)
risk    <- risk.model$risk.regular.w
c.index <- Hmisc::rcorr.cens(risk, y)[1] # TODO: implement this as function

# estimated ATE
risk.model$ate.hat - abs(ate) # accurate!

# observed and predicted beneefit by quartile group
benefits <- get.benefits(risk.model)
benefits$group.observed.benefit
benefits$group.predicted.benefit

# coefficients of 2nd stage
risk.model$coefficients.stage2

### 2. effect modeling ----
effect.model <- effect.modeling(x = x, w = w, y = y, alpha = 1) 

effect.model$ate.hat - abs(ate) # accurate!

calibration.plot(effect.model)

