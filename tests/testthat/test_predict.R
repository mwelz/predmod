# load package
library("predmod", quietly = TRUE)

set.seed(1)
n <- 200
p <- 5
w <- rbinom(n, 1, 0.5)
X <- matrix(runif(n * p), n)
colnames(X) <- paste0("var",1:p)
pi <- plogis(X[,1] + X[,2] + w * X[,2])
status <- rbinom(n, 1, pi)

newn <- 500
neww <- rbinom(newn, 1, 0.5)
newX <- matrix(runif(newn * p), newn)
newz <- runif(newn)
colnames(newX) <- paste0("var",1:p)

# baseline risk
BR <- baseline_risk(X, status)

# risk & effect model
RM <- risk_model(X, status, w)
EM <- effect_model(X, status, w)

test_that("check that predict methods works as intended",{
  
  # baseline risk
  b0 <- BR$risk
  b1 <- predict.baseline_risk(object = BR, newX = X, shrunk = FALSE)
  expect_equal(b0, b1)
  b2 <- predict.baseline_risk(object = BR, newX = newX, shrunk = FALSE) 
  expect_false(isTRUE(all.equal(b0, b2))) # different data => different results
  
  # risk model
  b1 <- predict.risk_model(RM, neww = w, newz = qlogis(RM$risk$baseline))[,"benefit_absolute"]
  b2 <- predict.risk_model(RM, neww = w, newz = RM$inputs$z)[,"benefit_absolute"]
  b0 <- RM$benefits$absolute
  expect_equal(b0, b1)
  expect_equal(b0, b2)
  b3 <- predict.risk_model(RM, neww = w, newz = RM$inputs$z)
  b4 <- predict.risk_model(RM, neww = neww, newz = newz, shrunk = FALSE) 
  expect_false(isTRUE(all.equal(b3, b4))) # different data => different results
  
  # effect model
  b1 <- predict.effect_model(EM, neww=w, newX=X)[,"benefit_absolute"]
  b0 <- EM$benefits$absolute
  expect_equal(b0, b1)
  b2 <- predict.effect_model(EM, neww=w, newX=X)
  b3 <- predict.effect_model(EM, neww=neww, newX=newX) 
  expect_false(isTRUE(all.equal(b2, b3))) # different data => different results
  
})