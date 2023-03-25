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
x <- causal_forest(X, status, w, num_trees = 100)
z <- runif(n)
GRF <- lapply(1:3, function(...) x)

test_that("check that get_benefits_causal_forest() behaves as expected",{
  
  expect_error(get_benefits_causal_forest(x), NA)
  expect_error(calibration_plot_causal_forest(x), NA)
  
  p0 <- get_benefits_causal_forest(x, baseline_risk = z)
  p1 <- get_benefits_causal_forest_imputation(GRF, baseline_risk = lapply(1:3, function(...) z))
  expect_equal(p0$absolute_predicted_benefit, p1$absolute_predicted_benefit)
  
  p0 <- calibration_plot_causal_forest(x, baseline_risk = z)
  p1 <- calibration_plot_causal_forest_imputation(GRF, baseline_risk = lapply(1:3, function(...) z))
  expect_equal(p0$data, p1$data)
  
})
