rm(list = ls()) ; cat("\014")

# load functions
source(paste0(getwd(), "/funs/linear-models/effect-modeling.R"))
source(paste0(getwd(), "/funs/plotmakers/plotmakers.R"))


### 0. Data generation ---- 
set.seed(22)
n <- 1000
p <- 5

# treatment assignment, 50% treatment, 50% control.
w <- rbinom(n, 1, 0.5) 

# covariates for outcome variable (y)
X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
colnames(X) <- paste0("nam.", 1:p)

# coefficients (including an intercept)
theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)

# compute Pr(Y = 1 | X) for each individual (with noise)
pi0 <- plogis(as.numeric(cbind(1,X) %*% theta))
pi1 <- plogis(as.numeric(cbind(1,X) %*% theta) + 1)

# create binary outcomes
y0 <- rbinom(n, 1, pi0)
y1 <- rbinom(n, 1, pi1)
y  <- ifelse(w == 1, y1, y0) # observed outcome


### 1. effect modeling ----
# arguments
alpha = 1
sig.level = 0.05
interacted.variables = colnames(X) 
retained.variables = c(colnames(X), "w.nam.1") # A string array of variable names in Z that shall always be retained (also works on interaction variables). If NULL, then no restriction applies. Note that treatment assignment w will always be retained by the function.
significance.level = 0.05
prediction.timeframe = NULL
lifeyears = NULL

# run effect model
em <- effect.modeling(X = X, y = y, w = w, 
                       interacted.variables = interacted.variables, 
                       alpha = alpha, lifeyears = lifeyears, 
                       prediction.timeframe = prediction.timeframe,
                       retained.variables = retained.variables, 
                       significance.level = significance.level)

em$effect.model$model.building$selecton.process$significance.tests_coefficients

calibration.plot(em)
subgroup.plot(em, X[,1])

em$effect.model$summary
em$average.treatment.effect
em$C.statistics
