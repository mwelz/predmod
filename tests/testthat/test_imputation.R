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
baseline_risk <- runif(n)
em_br <- effect_model(X, status, w, baseline_risk = baseline_risk)
em <- effect_model(X, status, w)
rm_br <- risk_model(X, status, w, z = runif(n))
rm <- risk_model(X, status, w)


m <- 3
BR <- lapply(1:m, function(...) br)
EM_br <- lapply(1:m, function(...) em_br)
EM <- lapply(1:m, function(...) em)
RM_br <- lapply(1:m, function(...) rm_br)
RM <- lapply(1:m, function(...) rm)

## further tests on newX etc. can be found in the tester for the non-imputation get_benefits() [if that one works, also the imputation generalization will work by construction]

test_that("check that the general imputation accounters do what they should",{
  
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
  


test_that("check that the benefits imputation accounters do what they should",{
  
  ## effect models 
  ord <- get_benefits(em)
  imp <- get_benefits_imputation(EM)
  expect_equal(ord$predicted_benefit, imp$predicted_benefit)
  expect_equal(ord$observed_benefit, imp$observed_benefit)
  
  ord <- get_benefits(em_br)
  imp <- get_benefits_imputation(EM_br)
  expect_equal(ord$predicted_benefit, imp$predicted_benefit)
  expect_equal(ord$observed_benefit, imp$observed_benefit)
  
  ord <- get_benefits(em, baseline_risk = baseline_risk)
  imp <- get_benefits_imputation(EM, baseline_risk = lapply(1:m, function(...) baseline_risk ))
  expect_equal(ord$predicted_benefit, imp$predicted_benefit)
  expect_equal(ord$observed_benefit, imp$observed_benefit)
  
  ## risk models
  ord <- get_benefits(rm)
  imp <- get_benefits_imputation(RM)
  expect_equal(ord$predicted_benefit, imp$predicted_benefit)
  expect_equal(ord$observed_benefit, imp$observed_benefit)
  
  ord <- get_benefits(rm_br, baseline_risk = baseline_risk)
  imp <- get_benefits_imputation(RM_br, baseline_risk = lapply(1:m, function(...) baseline_risk ))
  expect_equal(ord$predicted_benefit, imp$predicted_benefit)
  expect_equal(ord$observed_benefit, imp$observed_benefit)
  
  ord <- get_benefits(rm, baseline_risk = baseline_risk)
  imp <- get_benefits_imputation(RM, baseline_risk = lapply(1:m, function(...) baseline_risk ))
  expect_equal(ord$predicted_benefit, imp$predicted_benefit)
  expect_equal(ord$observed_benefit, imp$observed_benefit)
  
})


test_that("check that the calibration plot imputation accounters do what they should",{
  
  ## those are direct derivatives of the get benefits functions, so they should be fine
  
  ## effect models 
  ord <- calibration_plot(em)
  imp <- calibration_plot_imputation(EM)
  expect_equal(ord$data, imp$data)

  ord <- calibration_plot(em_br)
  imp <- calibration_plot_imputation(EM_br)
  expect_equal(ord$data, imp$data)
  
  ord <- calibration_plot(em, baseline_risk = baseline_risk)
  imp <- calibration_plot_imputation(EM, baseline_risk = lapply(1:m, function(...) baseline_risk ))
  expect_equal(ord$data, imp$data)
  
  ## risk models
  ord <- calibration_plot(rm)
  imp <- calibration_plot_imputation(RM)
  expect_equal(ord$data, imp$data)
  
  ord <- calibration_plot(rm_br, baseline_risk = baseline_risk)
  imp <- calibration_plot_imputation(RM_br, baseline_risk = lapply(1:m, function(...) baseline_risk ))
  expect_equal(ord$data, imp$data)
  
  ord <- calibration_plot(rm, baseline_risk = baseline_risk)
  imp <- calibration_plot_imputation(RM, baseline_risk = lapply(1:m, function(...) baseline_risk ))
  expect_equal(ord$data, imp$data)
  
})
