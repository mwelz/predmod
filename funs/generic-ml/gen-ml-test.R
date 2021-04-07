rm(list = ls()) ; cat("\014")

source(paste0(getwd(), "/funs/generic-ml/generic-ml-estimation-funs.R"))

logistic <- function(x) 1 / ( 1 + exp(-x))

set.seed(1)
num.obs  <- 500
num.vars <- 5

# treatment assignment (assume RCT)
D <- rbinom(num.obs, 1, 0.5) 

# covariates
Z <- mvtnorm::rmvnorm(num.obs, mean = rep(0, num.vars), sigma = diag(num.vars))

# coefficients (including an intercept)
theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)

# compute Pr(Y = 1 | X) for each individual (with noise)
pi0 <- logistic(as.numeric(cbind(1, Z) %*% theta))

# relative constant risk reduction of 30%
scaling <- 0.7
pi1 <- pi0 * scaling 
hte <- pi1 - pi0

# true ATE
ate <- mean(pi1) - mean(pi0)

# create binary outcomes
Y0 <- rbinom(num.obs, 1, pi0)
Y1 <- rbinom(num.obs, 1, pi1)
Y  <- ifelse(D == 1, Y1, Y0) # observed outcome

#####################

# require: D, Z, Y, the learners, significance level
# initialize
N     <- length(Y)
N.set <- 1:N

proportion.in.main.set = 0.5 # argument

### step 1: compute propensity scores ----
propensity.scores.obj <- propensity.score(Z = Z, D = D, learner = "glm")
propensity.scores     <- propensity.scores.obj$propensity.scores

### step 2: randomly split sample into main set and auxiliary set A
M.set <- sort(sample(x = N.set, size = floor(proportion.in.main.set * N), replace = FALSE),
              decreasing = FALSE)
A.set <- setdiff(N.set, M.set)



### step 2a: learn proxy predictors by using the auxiliary set

# get the proxy baseline estimator for the main sample
proxy.baseline.obj <- baseline.proxy.estimator(Z = Z, D = D, Y = Y, 
                                               auxiliary.sample = A.set, learner = "glm")
proxy.baseline     <- proxy.baseline.obj$baseline.predictions.main.sample

# get the proxy estimator of the CATE for the main sample
proxy.cate.obj <- CATE.proxy.estimator(Z = Z, D = D, Y = Y,
                                   auxiliary.sample = A.set, learner = "glm",
                                   proxy.baseline.estimates = proxy.baseline.obj$baseline.predictions.full.sample)
proxy.cate <- proxy.cate.obj$CATE.predictions.main.sample
