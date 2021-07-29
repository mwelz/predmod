rm(list = ls()) ; cat("\014")

source(paste0(getwd(), "/funs/grf/grf-funs.R"))
source(paste0(getwd(), "/funs/plotmakers/plotmakers.R"))
source(paste0(getwd(), "/funs/imputation/imputation.R"))


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
grf.model.imputed        <- rep(list(NA_real_), m)
names(grf.model.imputed) <- paste0("imputed.model_", 1:m)


for(i in 1:m){
  
  grf.model.imputed[[i]] <- 
    grf.modeling(X = imputed.datasets[[i]]$X, 
                    y = imputed.datasets[[i]]$y,
                    w = imputed.datasets[[i]]$w, num.trees = 500) 
  
} # FOR m

grf.obj <- grf.modeling_imputation.accounter(grf.model.imputed)

calibration.plot.grf_imputation.accounter(grf.model.imputed)
