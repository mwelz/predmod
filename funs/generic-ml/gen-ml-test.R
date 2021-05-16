rm(list = ls()) ; cat("\014")
library(glmnet) # apparently needs to be loaded (?) TODO
library(ggplot2)


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
learners.genericML <- "glm" #'mlr3::lrn("ranger", num.trees = 50)'# c('glm') #, 'tree', 'mlr3::lrn("ranger", num.trees = 50)') plot below is weird if you uncomment!
learner.propensity.score <- 'mlr3::lrn("glmnet", lambda = 0, alpha = 1)' # non-penalized logistic regression
num.splits <- 2
significance.level <- 0.05
store.splits <- FALSE
store.learners <- FALSE

# TODO: The TODOs in genericML()
genML <- genericML(Z = Z, D = D, Y = Y, 
                   learner.propensity.score = learner.propensity.score, 
                   learners.genericML = learners.genericML,
                   num.splits = num.splits,
                   Z.clan = Z.clan,
                   quantile.cutoffs = quantile.cutoffs,
                   proportion.in.main.set = proportion.in.main.set, 
                   significance.level = significance.level,
                   store.splits = store.splits,
                   store.learners = store.learners)


# analyze
genML$VEIN$best.learners$GATES # difference is insignificant, so no hetero
genML$VEIN$best.learners$BLP  # beta2 is insignificant, so no hetero
genML$VEIN$best.learners$CLAN$z1 # there seems to be hetero along z1



# plot
df <- data.frame(gates = genML$VEIN$best.learners$GATES[, "Estimate"],
                 ci.lower = genML$VEIN$best.learners$GATES[, "CI lower"],
                 ci.upper = genML$VEIN$best.learners$GATES[, "CI upper"],
                 group = 1:5)

ggplot(mapping = aes(x = group,
                     y = gates), data = df) +
  geom_hline(aes(yintercept = genML$VEIN$best.learners$BLP["beta.1", "Estimate"],
                 color = "ATE"),
             linetype = "dashed") +
  geom_hline(aes(yintercept = genML$VEIN$best.learners$BLP["beta.1", "CI lower"],
                 color = "ATE(90% CB)"),
             linetype = "dashed")  +
  geom_hline(yintercept = genML$VEIN$best.learners$BLP["beta.1", "CI upper"],
             linetype = "dashed", color = "red") +
  geom_point(aes(color = "GATES with 90% CB"), size = 3) +
  geom_errorbar(mapping = aes(ymin = ci.lower,
                              ymax = ci.upper)) +
  theme_light() +
  ylab("Treatment Effect") +
  xlab("Group by HTE Score") +
  scale_colour_manual(values = c("blue","red", "black")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom") + 
  guides(color = guide_legend(override.aes = list(
    linetype = 0, size = 4, shape = 15, alpha = 1))
  )
  
