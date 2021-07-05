rm(list = ls()) ; cat("\014")

# load functions
source(paste0(getwd(), "/funs/linear-models/risk-modeling.R"))
source(paste0(getwd(), "/funs/linear-models/effect-modeling.R"))
source(paste0(getwd(), "/funs/grf/grf-funs.R"))
source(paste0(getwd(), "/funs/plotmakers/plotmakers.R"))
source(paste0(getwd(), "/funs/imputation/imputation.R"))

# TODO: write documentation in risk modeling; effetc modeling etc
# TODO: subgroup plot imputation
# TODO: write function for input of effect, grf, risk modeling with imputation (see below for inspiration)
# TODO: grf imputer accounter (slightly different, as we have stderrs in returned objects already)


### 1. make data ----
set.seed(2)
n <- 1000
p <- 5

# treatment assignment, 50% treatment, 50% control.
w <- rbinom(n, 1, 0.5) 

# covariates for outcome variable (y)
X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
colnames(X) <- paste0("nam", 1:p)

# coefficients (including an intercept)
theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)

# compute Pr(Y = 1 | X) for each individual (with noise)
eps <-  rnorm(n, mean = 0, sd = 0.5)
pi0 <- plogis(as.numeric(cbind(1,X) %*% theta) + eps)

#assume a true relative constant risk reduction of 30%
scaling <- 0.7
pi1 <- pi0 * scaling 

# create binary outcomes
y0 <- rbinom(n, 1, pi0)
y1 <- rbinom(n, 1, pi1)
y  <- ifelse(w == 1, y1, y0) # observed outcome

# make some missing values in X
X[1:2, ] <- NA


### 2. impute data ----
m <- 10
imputed.datasets <- multiple.imputation(X = X, y = y, w = w, m = m, k = 5)


### 3. run model k times ----
# initialize
predictive.model.imputed        <- rep(list(NA_real_), m)
names(predictive.model.imputed) <- paste0("imputed.model_", 1:m)

for(i in 1:m){
  
  # fill in data in list
  alpha <- 1
  prediction.timeframe = NULL
  lifeyears = NULL
  
  predictive.model.imputed[[i]] <- 
    risk.modeling(X = imputed.datasets[[i]]$X, 
                  y = imputed.datasets[[i]]$y,
                  w = imputed.datasets[[i]]$w,
                  alpha = alpha) # TODO: effect modeling doesn't really work here for some reason!
  
} # FOR m


### 4. account for imputation uncertainty ----
rm.imp <- risk.modeling_imputation.accounter(predictive.model.imputed)

cutoffs = c(0.25, 0.5, 0.75)
significance.level = 0.05
rm.ben <- get.benefits_imputation.accounter(predictive.model.imputed, cutoffs, significance.level)

library(ggplot2)
calibration.plot_imputation.accounter(predictive.model.imputed)