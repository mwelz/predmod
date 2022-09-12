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
br <- baseline_risk(X, status)
em_br <- effect_model(X, status, w, baseline_risk = runif(n))
em <- effect_model(X, status, w)
rm_br <- risk_model(X, status, w, z = runif(n))
rm <- risk_model(X, status, w)


m <- 3
BR <- lapply(1:m, function(...) br)
EM_br <- lapply(1:m, function(...) em_br)
EM <- lapply(1:m, function(...) em)
RM_br <- lapply(1:m, function(...) rm_br)
RM <- lapply(1:m, function(...) rm)


test_that("check that the imputation accounters so what they should",{
  
  obj <- impaccount_baseline_risk(BR)
  expect_equal(obj$coefficients$full, br$coefficients$full)
  expect_equal(obj$coefficients$reduced, br$coefficients$reduced)
  
  obj <- impaccount_effect_model(EM_br)
  expect_equal(obj$coefficients$reduced, em_br$coefficients$reduced)
  
  obj <- impaccount_effect_model(EM)
  expect_equal(obj$coefficients$baseline$reduced, em$coefficients$baseline$reduced)
  expect_equal(obj$coefficients$reduced, em$coefficients$reduced)
  
  obj <- impaccount_risk_model(RM_br)
  expect_equal(obj$coefficients$stage2$accepted, rm_br$coefficients$stage2$accepted)

  obj <- impaccount_risk_model(RM)
  expect_equal(obj$coefficients$stage2$accepted, rm$coefficients$stage2$accepted)
  
})
  

