rm(list = ls()) ; cat("\014")

source(paste0(getwd(), "/funs/cox/cox-effect-modeling.R"))

logistic <- function(x) 1 / (1 + exp(-x))

### 0.1. Data generation ----
set.seed(2)
n <- 200
p <- 5

# treatment assignment (assume RCT)
w <- rbinom(n, 1, 0.5) 

# covariates
x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))

# coefficients (including an intercept)
theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)

# compute Pr(Y = 1 | X) for each individual (with noise)
eps <-  rnorm(n, mean = 0, sd = 0.5)
pi0 <- logistic(as.numeric(cbind(1, x) %*% theta) + eps)

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


### 0.2. Generating lifeyears
# Assume base 5 Lifeyears per person, assuming higher risk for disease is correlated with lower life-years and that those who die from disease lose 10% of their lifeyears left
Base_LY   <- 5
LY_coeffs <- c(0.5, -0.3, 0.7, -0.1, 0.4) #reversed signs of the coeffs for diseaseas placeholder
LY_intermediate <- (Base_LY + x  %*% LY_coeffs) * ifelse(y == 1, 0.9, 1) 
# set LY of 0 or less to 0 LY
LY <- ifelse(LY_intermediate <=0, 0, LY_intermediate) 

X = x
colnames(X) = paste0("nam", 1:p)
prediction.timeframe = 10
alpha = 1
lifeyears = LY
interacted.variables = colnames(X)
retained.variables = NULL
significance.level = 0.05

em <- cox.effect.modeling(X, y, w, interacted.variables = interacted.variables, alpha = alpha, lifeyears = lifeyears, prediction.timeframe = prediction.timeframe, retained.variables = retained.variables, significance.level = significance.level)
