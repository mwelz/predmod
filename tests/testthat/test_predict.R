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
  
  # risk model
  b1 <- predict.risk_model(RM, neww = w, newz = qlogis(RM$risk$baseline))[,"benefit_absolute"]
  b2 <- predict.risk_model(RM, neww = w, newz = RM$inputs$z)[,"benefit_absolute"]
  b0 <- RM$benefits$absolute
  expect_equal(b0, b1)
  expect_equal(b0, b2)
  
  # effect model
  b1 <- predict.effect_model(EM, neww=w, newX=X)[,"benefit_absolute"]
  b0 <- EM$benefits$absolute
  expect_equal(b0, b1)
  
})