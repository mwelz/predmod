rm(list = ls()) ; cat("\014")

# load functions
source(paste0(getwd(), "/funs2/estimation-funs/effect-modeling.R"))

### 0. Data generation ---- 
set.seed(2)
n <- 1000
p <- 5

# treatment assignment, 50% treatment, 50% control.
w <- rbinom(n, 1, 0.5) 

# covariates for outcome variable (y)
X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
colnames(X) <- paste0("nam", 1:p)

# coefficients (including an intercept)
theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)

# compute Pr(Y = 1 | X) for each individual (with noise)
eps <-  rnorm(n, mean = 0, sd = 0.5)
pi0 <- plogis(as.numeric(cbind(1,X) %*% theta) + eps)

#assume a true relative constant risk reduction of 30%
scaling <- 0.7
pi1 <- pi0 * scaling 

# create binary outcomes
y0 <- rbinom(n, 1, pi0)
y1 <- rbinom(n, 1, pi1)
y  <- ifelse(w == 1, y1, y0) # observed outcome


### 1. effect modeling ----
# arguments
alpha = 1
sig.level = 0.05
interacted.variables = colnames(X) 
retained.variables = c(colnames(X)) # A string array of variable names in Z that shall always be retained (also works on interaction variables). If NULL, then no restriction applies. Note that treatment assignment w will always be retained by the function.
significance.level = 0.05

# run effect model
em <- effect.modeling(X = X, y = y, w = w, 
                       interacted.variables = interacted.variables, 
                       alpha = alpha, 
                       retained.variables = retained.variables, 
                       significance.level = significance.level)

em$effect.model$summary
em$ate.hat
em$c.index.outcome
em$c.index.benefit

# TODO: incorporate lifeyears and predictiontimeframe as arguments
# TODO: make baseline risk return a probability (i.e. "response") instead of:
# baseline.mod <- risk.model.stage1(X = X, y = y, alpha = alpha)
# basline.risk <- transform.to.probability(baseline.mod$lp) 
# TODO: add return for baseline.model
# TODO: maybe function for the C benefits?