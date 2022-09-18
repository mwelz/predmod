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
br <- as.numeric(predict.baseline_risk(BR, newX = X, shrunk = FALSE))

# risk & effect model
RM <- risk_model(X = X, status = status, w = w, z = br)
EM <- effect_model(X, status, w, baseline_risk = br)


test_that("check that concordance functions work as intended for risk model",{
  
  C <- concordance(RM)
  
  # 'shrunk' is a dead argument in risk models
  pr      <- predict.risk_model(RM, neww = w, newz = br, shrunk = FALSE)
  risk    <- pr[,"risk_regular"]
  predben <- pr[,"benefit_absolute"]
  
  # calculate concordance
  cben <- C_benefit(y = status, w = w, pred_ben = predben)
  cout <- C_outcome(y = status, risk = risk)
  
  # test
  expect_equal(C$outcome, cout)
  expect_equal(C$benefit, cben)
})


test_that("check that concordance functions work as intended for effect model",{
  
  C <- concordance(EM)
  
  pr     <- predict.effect_model(EM, neww = w, newX = X, shrunk = FALSE)
  risk    <- pr[,"risk_regular"]
  predben <- pr[,"benefit_absolute"]
  
  # calculate concordance
  cben <- C_benefit(y = status, w = w, pred_ben = predben)
  cout_risk <- C_outcome(y = status, risk = risk)
  cout_base <- C_outcome(y = status, risk = br)
  
  
  # test
  expect_equal(C$outcome_baseline, cout_base)
  expect_equal(C$outcome, cout_risk)
  expect_equal(C$benefit, cben)
})
  