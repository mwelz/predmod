rm(list = ls()) ; cat("\014")

library(glmnet)
library(ggplot2)
source(paste0(getwd(), "/funs/estimation-funs.R"))


### 0 Data generation ----
set.seed(2)
n <- 10000
p <- 5

# ate parameter
theta <- -2

# treatment assignment (assume RCT)
w <- rbinom(n, 1, 0.5) 

# covariates
X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))

# coefficients (including an intercept)
beta <- c(1, 2, -3, 0, 0, 2)

# compute Pr(Y = 1 | X) for each individual (with noise)
#eps     <-  rnorm(n, mean = 0, sd = 0.5)
y0.star <- as.numeric(cbind(1, X) %*% beta) #+ eps
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
y0 <- rbinom(n, 1, pi0)
y1 <- rbinom(n, 1, pi1)
y  <- ifelse(w == 1, y1, y0) # observed outcome

#hist(pi0)
#hist(pi1)

# GLM according to DGP
fit.logistic <- glm(y ~., data = data.frame(y, X, w), family = binomial)
fit.logistic$coefficients
c(beta, theta) # good performance

# Use regularized regression Should kick out X3 and X4
fit.lasso    <- glmnet::cv.glmnet(cbind(X, w), y, family = "binomial", alpha = 1)
coefs.lasso  <- glmnet::coef.glmnet(fit.lasso, s = "lambda.min")
as.numeric(coefs.lasso)
c(beta, theta) # good performance

# now use retained variables in a normal non-penalized regression
kept.vars       <- coefs.lasso@i 
if(0 %in% kept.vars){
  kept.vars <- kept.vars[-which(kept.vars == 0)]
}
X.reduced <- data.frame(X, w)[,kept.vars]

fit.reduced <- glm(y ~., data = data.frame(y, X.reduced), family = binomial)
fit.reduced$coefficients
c(beta, theta) # good performance

# effect model
em <- effect.modeling(X, w, y, alpha = 1)
em$effect.model$coefficients
c(beta, theta) # good performance
calibration.plot(em, relative = TRUE)
calibration.plot(em, relative = FALSE)
