rm(list = ls()) ; cat("\014")

# load the helper functions
source(paste0(getwd(), "/tester-funs.R"))

# define the logistic function
logistic <- function(x) 1 / (1 + exp(-x))

### 0. Data generation ----
set.seed(2)
n <- 10000
p <- 5

# treatment assignment (assume RCT)
w <- rbinom(n, 1, 0.5) 

# true treatment effect
tau <- -1

# covariates
x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))

# coefficients (including an intercept)
theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)

# compute Pr(Y = 1 | X) for each individual (with noise)
eps <-  rnorm(n, mean = 0, sd = 0.5)
#pi1 <- logistic(tau + as.numeric(cbind(1, x) %*% theta) + eps)
pi1 <- logistic(as.numeric(cbind(1, x) %*% theta) + eps) * 0.4 # relatively 20% reduction. Absolutely it will differ.
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

ben <- get.benefits(risk.model)
ben$group.predicted.benefit

# calibration plot with the 4 quartiles
calibration.plot(risk.model)

# c index (based on second stage risk probabilities)
risk    <- risk.model$risk.regular.w
(c.index <- Hmisc::rcorr.cens(risk, y)[1]) # TODO: implement this as function

# estimated ATE
risk.model$ate.hat - abs(ate) # accurate!

# coefficients of 2nd stage
risk.model$coefficients.stage2

# plot the grouped 1st variable
subgroup.plot(risk.model, x = x[,1], quantile.nam = FALSE)

plot(x[,1], risk.model$predicted.benefit, ylab = "Predicted benefit")




### 2. effect modeling ----
effect.model <- effect.modeling(x = x, w = w, y = y, alpha = 1) 
Hmisc::rcorr.cens(effect.model$risk.regular.w, y)[1]
effect.model$ate.hat - abs(ate) # accurate!
calibration.plot(effect.model)

effect.model$effect.model$formula
round(effect.model$effect.model$summary, 4)
theta
tau

subgroup.plot(effect.model, x = x[,1], quantile.nam = FALSE)
plot(x[,1], effect.model$predicted.benefit, ylab = "Predicted benefit")




# TODO: summary statistics (count) in subgroup.plot, order of x-axis and groups
# TODO: make the groups of observed benefit positive!! (see calibration plot!)
# TODO: X vs x in arguments
