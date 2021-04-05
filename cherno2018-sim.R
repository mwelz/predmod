rm(list = ls()) ; cat("\014")
# based on https://docs.doubleml.org/r/stable/articles/DoubleML.html

# surpress messages from mlr3 package during fitting
lgr::get_logger("mlr3")$set_threshold("warn")

# initialize 
set.seed(31441)

n        <- 10000
p        <- 5
R        <- 100
coverage <- rep(NA_real_, R)

logistic <- function(x) 1 / ( 1 + exp(-x))

for(r in 1:R){
 
  
  # treatment assignment (assume RCT)
  w <- rbinom(n, 1, 0.5) 
  
  # covariates
  x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
  
  # coefficients (including an intercept)
  theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)
  
  # compute Pr(Y = 1 | X) for each individual (with noise)
  pi0 <- logistic(as.numeric(cbind(1, x) %*% theta))
  
  # relative constant risk reduction of 30%
  scaling <- 0.7
  pi1 <- pi0 * scaling 
  hte <- pi1 - pi0
  
  # true ATE
  ate <- mean(pi1) - mean(pi0)
  
  # create binary outcomes
  y0 <- rbinom(n, 1, pi0)
  y1 <- rbinom(n, 1, pi1)
  y  <- ifelse(w == 1, y1, y0) # observed outcome
  
  
  
  # matrix interface to DoubleMLData
  dml_data_sim <- DoubleML::double_ml_data_from_matrix(X=x, y=y, d=w)
  
  # for our simulated data from a sparse partially linear model we use a Lasso regression model.
  ml_g_sim <- mlr3::lrn("regr.glmnet", lambda = sqrt(log(p)/(n)))
  ml_m_sim <-  mlr3::lrn("classif.glmnet", lambda = sqrt(log(p)/(n)))
  
  # estimate and fit model
  obj_dml_plr_sim <- DoubleML::DoubleMLIRM$new(dml_data_sim, ml_g=ml_g_sim, ml_m=ml_m_sim)
  obj_dml_plr_sim$fit()
  
  
  # evaluate model
  ci.lo <- obj_dml_plr_sim$confint(level = 0.95)[1]
  ci.up <- obj_dml_plr_sim$confint(level = 0.95)[2]

  
  coverage[r] <- (ci.lo <= ate) & (ate <= ci.up)
  
}

mean(coverage) # 0.96. Expected 95%, so that's okay

