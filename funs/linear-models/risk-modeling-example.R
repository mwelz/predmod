rm(list = ls()) ; cat("\014")

# load functions
source(paste0(getwd(), "/funs2/linear-models/risk-modeling.R"))
source(paste0(getwd(), "/funs2/plotmakers/plotmakers.R"))


# make data
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


# test the risk model
alpha <- 1
offset.linear.predictor <- TRUE
prediction.timeframe = NULL
lifeyears = NULL

RM <- risk.modeling(X = X, y = y, w = w, alpha = alpha, 
                    offset.linear.predictor = offset.linear.predictor, 
                    lifeyears = lifeyears, prediction.timeframe = prediction.timeframe)

calibration.plot(RM)
subgroup.plot(RM, X[,1])

RM$models$coefficients.stage2
RM$average.treatment.effect
RM$C.statistics
