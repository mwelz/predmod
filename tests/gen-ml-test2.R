rm(list = ls()) ; cat("\014")
library(glmnet) # apparently needs to be loaded (?) TODO
library(ggplot2)


source(paste0(getwd(), "/funs/generic-ml/generic-ml-estimation-funs.R"))
source(paste0(getwd(), "/funs/generic-ml/generic-ml-auxiliary-funs.R"))


logistic <- function(x) 1 / ( 1 + exp(-x))

set.seed(1)
num.obs  <- 500
num.vars <- 5

# ate parameter
theta <- -2

# treatment assignment (assume RCT)
D <- rbinom(num.obs, 1, 0.5) 

# covariates
Z <- mvtnorm::rmvnorm(num.obs, mean = rep(0, num.vars), sigma = diag(num.vars))
colnames(Z) <- paste0("z", 1:num.vars)

# coefficients (including an intercept)
beta <- c(1, 2, -3, 0, 0, 2)

# compute Pr(Y = 1 | X) for each individual (with noise)
#eps     <-  rnorm(n, mean = 0, sd = 0.5)
y0.star <- as.numeric(cbind(1, Z) %*% beta) #+ eps
y1.star <- theta + y0.star

# experimental: use standardized logit link
mu <- 0 # mean(y0.star)
s <- 1 #sqrt(var(y0.star))
pi0 <- plogis(y0.star, location = mu, scale = s)
pi1 <- plogis(y1.star, location = mu, scale = s)

# hte
hte <- pi1 - pi0

# true ATE
ate <- mean(pi1) - mean(pi0)

# create binary outcomes
y0 <- rbinom(num.vars, 1, pi0)
y1 <- rbinom(num.vars, 1, pi1)
Y  <- ifelse(D == 1, y1, y0) # observed outcome


#######################

# arguments: 
quantile.cutoffs         <- c(0.2, 0.4, 0.6, 0.8) # for the GATES grouping of S (argument)
proportion.in.main.set   <- 0.5 # argument
Z.clan                   <- NULL # argument. The matrix of variables that shall be considered in CLAN
learners.genericML       <- c("glm", 'mlr3::lrn("ranger", num.trees = 100)') 
learner.propensity.score <- 'mlr3::lrn("glmnet", lambda = 0, alpha = 1)' # non-penalized logistic regression
num.splits               <- 100
significance.level       <- 0.05
store.splits             <- FALSE
store.learners           <- FALSE


### step 1: compute propensity scores ----
propensity.scores.obj <- propensity.score(Z = Z, D = D, 
                                          learner = make.mlr3.string(learner.propensity.score, 
                                                                     regr = FALSE))
propensity.scores     <- propensity.scores.obj$propensity.scores

### step 2: for each ML method, do the generic ML analysis ----
learners = learners.genericML

# initialize
generic.targets <- initializer.for.splits(Z = Z, Z.clan = Z.clan, 
                                          learners = learners, num.splits = num.splits, 
                                          quantile.cutoffs = quantile.cutoffs)

num.vars.in.Z.clan <- ifelse(is.null(Z.clan), ncol(Z), ncol(Z.clan))
genericML.by.split <- list()
N     <- length(Y)
N.set <- 1:N

if(store.splits) splits.mat <- matrix(NA_character_, N, num.splits)

s=1
i=1
# loop over the sample splits
for(s in 1:num.splits){
  
  # perform sample splitting into main set and auxiliary set
  M.set <- sort(sample(x = N.set, size = floor(proportion.in.main.set * N), replace = FALSE),
                decreasing = FALSE)
  A.set <- setdiff(N.set, M.set)
  
  if(store.splits){
    
    splits.mat[M.set, s] <- "M"
    splits.mat[A.set, s] <- "A"
    
  } # IF
  
  
  # loop over the learners
  for(i in 1:length(learners)){
    
    learner = learners[i]
    
   
    ### step 1: input checks ---- 
    if(is.null(Z.clan)) Z.clan <- Z # if no input provided, set it equal to Z
    
    ### step 2a: learn proxy predictors by using the auxiliary set ----
    
    # get the proxy baseline estimator for the main sample
    proxy.baseline.obj <- baseline.proxy.estimator(Z = Z, D = D, Y = Y, 
                                                   auxiliary.sample = A.set, 
                                                   learner = make.mlr3.string(learner, regr = TRUE))
    proxy.baseline     <- proxy.baseline.obj$baseline.predictions.main.sample
    
    # get the proxy estimator of the CATE for the main sample
    proxy.cate.obj <- 
      CATE.proxy.estimator(Z = Z, D = D, Y = Y,
                           auxiliary.sample = A.set, 
                           learner = make.mlr3.string(learner, regr = TRUE),
                           proxy.baseline.estimates = proxy.baseline.obj$baseline.predictions.full.sample)
    proxy.cate <- proxy.cate.obj$CATE.predictions.main.sample
    
    
    ### step 2b: estimate BLP parameters by OLS (TODO: HT transformation!) ----
    blp.obj <- get.BLP.params.classic(D = D[M.set], Y = Y[M.set],
                                      propensity.scores = propensity.scores[M.set],
                                      proxy.baseline = proxy.baseline, 
                                      proxy.cate = proxy.cate, 
                                      significance.level = significance.level)
    
    
    # prepare weights
    weights <- 1 / (propensity.scores[M.set] * (1 - propensity.scores[M.set]))
    
    # prepare covariate matrix
    X <- data.frame(B = proxy.baseline, 
                    S = proxy.cate,
                    beta.1 = D - propensity.scores[M.set], 
                    beta.2 = (D - propensity.scores[M.set]) * (proxy.cate - mean(proxy.cate))) 
    
    # fit weighted linear regression by OLS
    blp.obj <- lm(Y ~., data = data.frame(Y=Y[M.set], X[M.set,]), weights = weights)
    
    # extract coefficients
    coefficients     <- summary(blp.obj)$coefficients
    
    # inference on beta2: test the null that it is 1) = 0, 2) = 1.
    beta2.inference <- matrix(NA_real_, 2, 2)
    rownames(beta2.inference) <- c("H0: beta.2 = 0", "H0: beta.2 = 1")
    colnames(beta2.inference) <- c("t value", "Pr(>|t|)")
    beta2.inference[1,] <- coefficients["beta.2", c("t value", "Pr(>|t|)")]
    beta2.inference[2, "t value"] <- 
      (coefficients["beta.2", "Estimate"] - 1) / coefficients["beta.2", "Std. Error"]  
    beta2.inference[2, "Pr(>|t|)"] <- 
      2 * pt(beta2.inference[2, "t value"], df = blp.obj$df.residual, lower.tail = FALSE)
    
    # generic targets
    generic.targets <- coefficients[c("beta.1", "beta.2"), ]
    colnames(generic.targets) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    generic.targets[,"Pr(>|z|)"] <- 2 * pnorm(abs(generic.targets[,"z value"]), lower.tail = FALSE)
    ci.lo <- generic.targets[,"Estimate"] - qnorm(1-significance.level/2) * generic.targets[,"Std. Error"]
    ci.up <- generic.targets[,"Estimate"] + qnorm(1-significance.level/2) * generic.targets[,"Std. Error"]
    generic.targets <- cbind(generic.targets, ci.lo, ci.up)
    colnames(generic.targets) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "CI lower", "CI upper")
    generic.targets <- generic.targets[,c("Estimate", "CI lower", "CI upper", 
                                          "Std. Error", "z value", "Pr(>|z|)")]
    
    
    
    
    ### step 2c: estimate GATES parameters by OLS (TODO: HT transformation!) ----
    # group the proxy estimators for the CATE in the main sample by quantiles. TODO: intervals need to be [) instead of (]
    group.membership.main.sample <- quantile.group(proxy.cate, 
                                                   cutoffs = quantile.cutoffs, 
                                                   quantile.nam = TRUE) 
    
    gates.obj <- get.GATES.params.classic(D = D[M.set], Y = Y[M.set],
                                          propensity.scores = propensity.scores[M.set],
                                          group.membership.main.sample = group.membership.main.sample, 
                                          proxy.baseline = proxy.baseline, proxy.cate = proxy.cate,
                                          significance.level = significance.level)
    
    
    ### step 2d: estimate CLAN parameters in the main sample
    clan.obj <- get.CLAN.parameters(Z.clan.main.sample = Z.clan[M.set,], 
                                    group.membership.main.sample = group.membership.main.sample)
    
    
    ### step 2e: get parameters over which we maximize to find the "best" ML method ----
    best.obj <- best.ml.method.parameters(BLP.obj = blp.obj, 
                                          GATES.obj = gates.obj, 
                                          proxy.cate.main.sample = proxy.cate,
                                          group.membership.main.sample = group.membership.main.sample)
    
    
    
    
    
    generic.targets[[i]]$BLP[,,s]   <- generic.ml.obj$BLP$generic.targets
    generic.targets[[i]]$GATES[,,s] <- generic.ml.obj$GATES$generic.targets
    generic.targets[[i]]$best[,,s]  <- c(generic.ml.obj$best$lambda, generic.ml.obj$best$lambda.bar)
    
    if(store.learners){
      
      genericML.by.split[[learners[i]]][[s]] <- generic.ml.obj
      
    }
    
    for(j in 1:num.vars.in.Z.clan){
      generic.targets[[i]]$CLAN[[j]][,,s] <- generic.ml.obj$CLAN$generic.targets[[j]]
    }
    
  } # FOR learners
} # FOR num.splits

if(!store.learners) genericML.by.split <- NULL
if(!store.splits)   splits.mat <- NULL

