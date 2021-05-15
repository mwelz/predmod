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
colnames(Z) <- paste0("z", 1:num.vars)

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


#######################

# arguments: 
quantile.cutoffs       = c(0.25, 0.5, 0.75) # for the GATES grouping of S (argument)
proportion.in.main.set = 0.5 # argument
Z.clan                 = NULL # argument. The matrix of variables that shall be considered in CLAN
learners.genericML <- c('glm', 'tree', 'mlr3::lrn("ranger", num.trees = 50)')
learner.propensity.score <- "glm"
num.splits <- 4
significance.level <- 0.05

# TODO: The TODOs in genericML()
genML <- genericML(Z = Z, D = D, Y = Y, 
                   learner.propensity.score = learner.propensity.score, 
                   learners.genericML = learners.genericML,
                   num.splits = num.splits,
                   Z.clan = Z.clan,
                   quantile.cutoffs = quantile.cutoffs,
                   proportion.in.main.set = proportion.in.main.set, 
                   significance.level = significance.level)


# analyze
genML$VEIN$best.learners$GATES # difference is insignificant, so no hetero
genML$VEIN$best.learners$BLP  # beta2 is insignificant, so no hetero
genML$VEIN$best.learners$CLAN$z1 # there seems to be hetero along z1
