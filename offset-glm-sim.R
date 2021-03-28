rm(list = ls()) ; cat("\014")

set.seed(1)

p <- 2; n <- 10000
beta0 <- 1
beta1 <- c(-1,1)
beta <- c(beta0, beta1)
x <- matrix(rnorm(n*p), nrow = n)
eta <- beta0 + beta1[1] * x[,1] + beta1[2] * x[,2] # the linear predictor
c <- runif(n, min = -2, max = -1) # the offset

### logistic regression ----
logistic <- function(x) 1 / (1 + exp(-x)) # the logistic function (inverse of the logit function)

# y is binary
probs.nooffset <- logistic(eta) # because in GLM, eta = logit(probs) [logit-link] => probs = logistic(eta), as logistic is inverse of logit.
probs.offset <- logistic(eta + c) # the offset is added to the linear predictor (this is how GLM assumes it to be)

y.nooffset <- rbinom(n, 1, probs.nooffset)
y.offset <- rbinom(n, 1, probs.offset)

sum(y.nooffset) # observe: the offset causes the number of successes to decrease dramatically!
sum(y.offset)

data <- data.frame(y.nooffset, y.offset, x1 = x[,1], x2 = x[,2], c = c)

# non-offsetted logit on non-offsetted data
logit.model1 <- glm(y.nooffset ~ x1 + x2, data = data, family = "binomial")
logit.model1$coefficients - beta  # good performance

# non-offsetted logit on offsetted data: we expect a bias in the estimate of beta0!
logit.model2 <- glm(y.offset ~ x1 + x2, data = data, family = "binomial")
logit.model2$coefficients - beta  # as expected, the estimate for the intercept is biased

# offsetted logit on offsetted data: no bias anymore expected!
logit.model3 <- glm(y.offset ~ x1 + x2, data = data, family = "binomial", offset = c)
logit.model3$coefficients - beta # as expected, the performance is good again!

# repeat this exercise by using glmnet instead of glm:
logit.model.net.1 <- glmnet::glmnet(x = cbind(x[,1], x[,2]), y = y.nooffset, family = "binomial", 
                                    lambda = 0, alpha = 1)
as.numeric(coefficients(logit.model.net.1)) - beta # good (as expected)


logit.model.net.2 <- glmnet::glmnet(x = cbind(x[,1], x[,2]), y = y.offset, family = "binomial", 
                                    lambda = 0, alpha = 1)
as.numeric(coefficients(logit.model.net.2)) - beta # bad (as expected)


logit.model.net.3 <- glmnet::glmnet(x = cbind(x[,1], x[,2]), y = y.offset, family = "binomial", 
                                    lambda = 0, alpha = 1, offset = c)
as.numeric(coefficients(logit.model.net.3)) - beta # good (as expected)

### poisson regression ----
# keep in mind that y is a counting process now
mu.nooffset <- exp(eta) # because in GLM, eta = log(mu) [log-link] => mu = exp(eta)
mu.offset <- exp(eta + c)  # the offset is added to the linear predictor
y.nooffset <- rpois(n = n, lambda = mu.nooffset)
y.offset <- rpois(n = n, lambda = mu.offset)
data <- data.frame(y.nooffset, y.offset, x1 = x[,1], x2 = x[,2], c = c)

# non-offsetted poisson on non-offsetted data 
poi.model1 <- glm(y.nooffset ~ x1 + x2, data = data, family = "poisson")
poi.model1$coefficients - beta # good performance

# non-offsetted on offsetted data: bias in estimate of intercept expected
poi.model2 <- glm(y.offset ~ x1 + x2, data = data, family = "poisson")
poi.model2$coefficients - beta # as expected, the estimate for the intercept is biased

# offsetted logit on offsetted data: no bias anymore expected!
poi.model3 <- glm(y.offset ~ x1 + x2, data = data, family = "poisson", offset = c)
poi.model3$coefficients - beta # as expected, the performance is good again!

# repeat this exercise by using glmnet instead of glm:
poi.model.net.1 <- glmnet::glmnet(x = cbind(x[,1], x[,2]), y = y.nooffset, family = "poisson", 
                                    lambda = 0, alpha = 1)
as.numeric(coefficients(poi.model.net.1)) - beta # good (as expected)


poi.model.net.2 <- glmnet::glmnet(x = cbind(x[,1], x[,2]), y = y.offset, family = "poisson", 
                                  lambda = 0, alpha = 1)
as.numeric(coefficients(poi.model.net.2)) - beta # bad (as expected)


poi.model.net.3 <- glmnet::glmnet(x = cbind(x[,1], x[,2]), y = y.offset, family = "poisson", 
                                  lambda = 0, alpha = 1, offset = c)
as.numeric(coefficients(poi.model.net.3)) - beta # good (as expected)
