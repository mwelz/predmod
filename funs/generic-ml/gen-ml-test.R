rm(list = ls()) ; cat("\014")
library(glmnet) # apparently needs to be loaded (?) TODO

source(paste0(getwd(), "/funs/generic-ml/generic-ml-estimation-funs.R"))
source(paste0(getwd(), "/funs/generic-ml/generic-ml-auxiliary-funs.R"))


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

quantile.cutoffs <- c(0.25, 0.5, 0.75) # for the GATES grouping of S (argument)
proportion.in.main.set = 0.5 # argument
Z.clan <- NULL # argument. The matrix of variables that shall be considered in CLAN

### step 0: input checks ---- 
if(is.null(Z.clan)) Z.clan <- Z # if no inputprovided, set it equal to Z

### step 1: compute propensity scores ----
propensity.scores.obj <- propensity.score(Z = Z, D = D, learner = "glm")
propensity.scores     <- propensity.scores.obj$propensity.scores

### step 2: randomly split sample into main set and auxiliary set A ----
M.set <- sort(sample(x = N.set, size = floor(proportion.in.main.set * N), replace = FALSE),
              decreasing = FALSE)
A.set <- setdiff(N.set, M.set)


### step 2a: learn proxy predictors by using the auxiliary set ----

# get the proxy baseline estimator for the main sample
proxy.baseline.obj <- baseline.proxy.estimator(Z = Z, D = D, Y = Y, 
                                               auxiliary.sample = A.set, learner = "glm")
proxy.baseline     <- proxy.baseline.obj$baseline.predictions.main.sample

# get the proxy estimator of the CATE for the main sample
proxy.cate.obj <- 
  CATE.proxy.estimator(Z = Z, D = D, Y = Y,
                       auxiliary.sample = A.set, learner = "glm",
                       proxy.baseline.estimates = proxy.baseline.obj$baseline.predictions.full.sample)
proxy.cate <- proxy.cate.obj$CATE.predictions.main.sample


### step 2b: estimate BLP parameters by OLS (TODO: HT transformation!) ----
blp.obj <- get.BLP.params.classic(D = D[M.set], Y = Y[M.set],
                                  propensity.scores = propensity.scores.obj$propensity.scores[M.set])


### step 2c: estimate GATES parameters by OLS (TODO: HT transformation!) ----
# group the proxy estimators for the CATE in the main sample by quantiles. TODO: intervals need to be [) instead of (]
group.membership.main.sample <- quantile.group(proxy.cate, 
                                cutoffs = quantile.cutoffs, 
                                quantile.nam = TRUE) 

gates.obj <- get.GATES.params.classic(D = D[M.set], Y = Y[M.set],
                                      propensity.scores = propensity.scores.obj$propensity.scores[M.set],
                                      group.membership.main.sample = group.membership.main.sample)


### step 2d: estimate CLAN parameters in the main sample
clan.obj <- get.CLAN.parameters(Z.clan.main.sample = Z.clan[M.set,], 
                                group.membership.main.sample = group.membership.main.sample)

# comments on CLAN: If there are categorical variables, apply one-hot-encoding to Z.clan. The interpretation then becomes: Is there a factor that is overproportionally present in the least or most affected group?


### step 2e: get parameters over which we maximize to find the "best" ML method ----
best.obj <- best.ml.method.parameters(BLP.obj = blp.obj, 
                                      GATES.obj = gates.obj, 
                                      proxy.cate.main.sample = proxy.cate,
                                      group.membership.main.sample = group.membership.main.sample)

