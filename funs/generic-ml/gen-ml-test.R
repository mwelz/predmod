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

#####################

# require: D, Z, Y, the learners, significance level
# arguments: 
quantile.cutoffs       = c(0.25, 0.5, 0.75) # for the GATES grouping of S (argument)
proportion.in.main.set = 0.5 # argument
Z.clan                 = NULL # argument. The matrix of variables that shall be considered in CLAN
learners <- c('glm', 'tree')
num.splits <- 2


### step 1: compute propensity scores ----
propensity.scores.obj <- propensity.score(Z = Z, D = D, learner = "glm")
propensity.scores     <- propensity.scores.obj$propensity.scores

### step 2: for each ML method, do the generic ML analysis
## put the following in a function later

# initialize
generic.targets <- initializer.for.splits(Z = Z, Z.clan = Z.clan, 
                                          learners = learners, num.splits = num.splits, 
                                          quantile.cutoffs = quantile.cutoffs)

num.vars.in.Z.clan <- ifelse(is.null(Z.clan), ncol(Z), ncol(Z.clan))

for(s in 1:num.splits){
  for(i in 1:length(learners)){
    
    generic.ml.obj <- 
      get.generic.ml.for.given.learner(Z = Z, D = D, Y = Y, 
                                       propensity.scores = propensity.scores, 
                                       learner = learners[i], 
                                       Z.clan = Z.clan, 
                                       proportion.in.main.set = proportion.in.main.set, 
                                       quantile.cutoffs = quantile.cutoffs)
    
    generic.targets[[i]]$BLP[,,s]   <- generic.ml.obj$BLP$generic.targets
    generic.targets[[i]]$GATES[,,s] <- generic.ml.obj$GATES$generic.targets
    generic.targets[[i]]$best[,,s]  <- c(generic.ml.obj$best$lambda, generic.ml.obj$best$lambda.bar)
    
    for(j in 1:num.vars.in.Z.clan){
      generic.targets[[i]]$CLAN[[j]] <- generic.ml.obj$CLAN$generic.targets[[j]]
    }
    
  } # FOR learners
} # FOR num.splits
