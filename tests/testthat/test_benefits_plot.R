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
baseline_risk <- as.numeric(BR$risk)

# risk model (silently fits a duplicate of BR)
RM <- risk_model(X = X, status = status, w = w)

# effect model
EM <- effect_model(X = X, status = status, w = w)

# new data
newX <- matrix(runif(300 * p), 300) # new X
newz <- runif(300)
neww <- rbinom(300, 1, 0.5)
newstatus <- rbinom(300, 1, 0.5)

newX2 <- matrix(runif(n * p), n) # new X
newz2 <- runif(n)
neww2 <- rbinom(n, 1, 0.5)
newstatus2 <- rbinom(n, 1, 0.5)



test_that("check that get_benefits() behaves as expected when new data are passed",{
  
  # incorrect input dimensions
  expect_error(get_benefits(RM, neww = neww, newz = z))
  expect_error(get_benefits(RM, neww = w, newz = newz)) 
  expect_error(get_benefits(EM, neww = neww, newX = X))
  expect_error(get_benefits(EM, neww = w, newz = newX))
  
  # newX is redundant in RM and therefore shouldn't matter
  expect_equal(unlist(get_benefits(RM)),
               unlist(get_benefits(RM, newX = newX)))
  
  # newz is redundant in EM and therefore shouldn't matter
  expect_equal(unlist(get_benefits(EM)),
               unlist(get_benefits(EM, newz = newz)))
  
  # passing a baseline risk with different length should throw an error
  expect_error(get_benefits(RM, baseline_risk = newz))
  expect_error(get_benefits(EM, baseline_risk = newz))
  expect_error(get_benefits(RM, newstatus = status, neww = w, 
                              newz = BR$linear_predictor, baseline_risk = newz))
  expect_error(get_benefits(EM, newstatus = status, neww = w, 
                              newX = X, baseline_risk = newz))
  
})


test_that("check that benefits and plot are the same on default as when passing in-sample data",{
  
  expect_equal(unlist(get_benefits(EM)), 
               unlist(get_benefits(EM, newstatus = status, newX = X, neww = w)))
  
  expect_equal(unlist(get_benefits(RM)), 
               unlist(get_benefits(RM, newstatus = status, neww = w, 
                                   newz = BR$linear_predictor, baseline_risk = baseline_risk)))
  
  # different data => different results
  b0 <- get_benefits(EM)
  b1 <- get_benefits(EM, newX = newX2, newstatus = newstatus2, neww = neww2)
  expect_false(isTRUE(all.equal(b0, b1))) 
  b0 <- get_benefits(RM)
  b1 <- get_benefits(RM, newz = newz2, newstatus = newstatus2, neww = neww2)
  expect_false(isTRUE(all.equal(b0, b1))) 
  
  
  # do the same for the calibration plots, which are derivatives of get_benefits()
  p0 <- calibration_plot(EM)
  p1 <- calibration_plot(EM, newstatus = status, newX = X, neww = w)
  expect_equal(unlist(p0$data), 
               unlist(p1$data))
  
  p0 <- calibration_plot(RM)
  p1 <- calibration_plot(RM, newstatus = status, neww = w, 
                         newz = BR$linear_predictor, baseline_risk = baseline_risk)
  expect_equal(unlist(p0$data), 
               unlist(p1$data))

})
