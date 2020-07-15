rm(list = ls()) ; cat("\014")

# load the helper functions
source(paste0(getwd(), "/risk-modeling/rm-funs/rm-funs.R"))

# define the logistic function
logistic <- function(x) 1 / (1 + exp(-x))

### 0. Data generation
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

## predictive modeling approach:
alpha <- 1
colnames(x) <- paste0("z_", 1:ncol(x))
interaction.vars <- c("z_2", "z_3", "z_4")


foo <- effect.modeling(x, w, y, alpha, interaction.vars = interaction.vars)

# TODO: bugfix in risk modeling functions: the kept.vars index at 0 refers to the coefficient, so we do not have zero indexing! (like the way we did it here)
# standardization of X before Lasso; see ESL p. 63
