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
subset <- sample.int(n, floor(n/2))

# risk model
z <- runif(n)
set.seed(1)
RM <- risk_model(X = X, status = status, w = w, z = z)

# effect model
EM <- effect_model(X = X, status = status, w = w)

# new data
newX <- matrix(runif(100 * p), 100) # new X
newz <- runif(100)
neww <- rbinom(100, 1, 0.5)

test_that("check that ATE expects as behaved when new data are passed",{
  
  # incorrect input dimensions
  expect_error(average_treatment_effect(RM, neww = neww, newz = z))
  expect_error(average_treatment_effect(RM, neww = w, newz = newz)) 
  expect_error(average_treatment_effect(EM, neww = neww, newX = X))
  expect_error(average_treatment_effect(EM, neww = w, newz = newX))
  
  # newX is redundant in RM and therefore shouldn't matter
  expect_equal(average_treatment_effect(RM),
               average_treatment_effect(RM, newX = newX))
  
  # newz is redundant in EM and therefore shouldn't matter
  expect_equal(average_treatment_effect(EM),
               average_treatment_effect(EM, newz = newz))
  
})



test_that("check that ATE in RM is the same on default as when passing in-sample data",{
  

  expect_equal(average_treatment_effect(RM), 
               average_treatment_effect(RM, neww = w, newz = z))
  
  expect_equal(unname(average_treatment_effect(RM)[1]), 
               mean(RM$benefits$absolute))
  
  expect_equal(average_treatment_effect(RM, relative = TRUE), 
               average_treatment_effect(RM, neww = w, newz = z, relative = TRUE))
  
  expect_equal(unname(average_treatment_effect(RM, relative = TRUE)[1]), 
               mean(RM$benefits$relative))
  
  expect_equal(average_treatment_effect(RM, subset = subset), 
               average_treatment_effect(RM, neww = w, newz = z, subset = subset))
  
  expect_equal(unname(average_treatment_effect(RM, subset = subset)[1]), 
               mean(RM$benefits$absolute[subset]))
  
  expect_equal(average_treatment_effect(RM, subset = subset, relative = TRUE), 
               average_treatment_effect(RM, neww = w, newz = z, subset = subset, relative = TRUE))
  
  expect_equal(unname(average_treatment_effect(RM, subset = subset, relative = TRUE)[1]), 
               mean(RM$benefits$relative[subset]))
  
}) # TEST


test_that("check that ATE in EM is the same on default as when passing in-sample data",{
  
  expect_equal(average_treatment_effect(EM), 
               average_treatment_effect(EM, neww = w, newX = X))
  
  expect_equal(unname(average_treatment_effect(EM)[1]), 
               mean(EM$benefits$absolute))
  
  expect_equal(average_treatment_effect(EM, relative = TRUE), 
               average_treatment_effect(EM, neww = w, newX  = X, relative = TRUE))
  
  expect_equal(unname(average_treatment_effect(EM, relative = TRUE)[1]), 
               mean(EM$benefits$relative))
  
  expect_equal(average_treatment_effect(EM, subset = subset), 
               average_treatment_effect(EM, neww = w, newX = X, subset = subset))
  
  expect_equal(unname(average_treatment_effect(EM, subset = subset)[1]), 
               mean(EM$benefits$absolute[subset]))
  
  expect_equal(average_treatment_effect(EM, subset = subset, relative = TRUE), 
               average_treatment_effect(EM, neww = w, newX = X, subset = subset, relative = TRUE))
  
  expect_equal(unname(average_treatment_effect(EM, subset = subset, relative = TRUE)[1]), 
               mean(EM$benefits$relative[subset]))
  
}) # TEST
