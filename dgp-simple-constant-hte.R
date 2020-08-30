rm(list = ls()) ; cat("\014")

# load the helper functions
source(paste0(getwd(), "/risk-modeling/rm-funs/rm-funs.R"))

# define the logistic function
logistic <- function(x) 1 / (1 + exp(-x))

### 0. Data generation ----
set.seed(2)
n <- 10000
p <- 5

# treatment assignment (assume RCT)
w <- rbinom(n, 1, 0.5) 

# covariates
x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))

# coefficients (including an intercept)
theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)

# compute Pr(Y = 1 | X) for each individual (with noise)
eps <-  rnorm(n, mean = 0, sd = 0.5)
pi0 <- logistic(as.numeric(cbind(1, x) %*% theta) + eps)

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


### 1. risk modeling ----
risk.model <- risk.modeling(X = x, w = w, y = y, alpha = 1, offset.lp = TRUE)


pdf(file = paste0(getwd(), "/tester-plots/const-rm-calibration-relative.pdf"))
calibration.plot(risk.model, relative = TRUE, title = "Risk Model, Calibration: Predicted Relative Benefit")
dev.off()

pdf(file = paste0(getwd(), "/tester-plots/const-rm-calibration-absolute.pdf"))
calibration.plot(risk.model, relative = FALSE, title = "Risk Model, Calibration: Predicted Absolute Benefit")
dev.off()

pdf(file = paste0(getwd(), "/tester-plots/const-rm-subgroup-absolute.pdf"))
subgroup.plot(risk.model, x[,1], relative = FALSE)
dev.off()

pdf(file = paste0(getwd(), "/tester-plots/const-rm-subgroup-relative.pdf"))
subgroup.plot(risk.model, x[,1], relative = TRUE)
dev.off()


### 2. effect modeling ----
effect.model <- effect.modeling(x = x, w = w, y = y, alpha = 1) 
# 
# calibration.plot(effect.model, relative = TRUE, title = "Effect Model, Calibration: Predicted Relative Benefit")
# calibration.plot(effect.model, relative = FALSE, title = "Effect Model, Calibration: Predicted Absolute Benefit")
# 
# subgroup.plot(effect.model, x[,1], relative = TRUE)
# subgroup.plot(effect.model, x[,1], relative = FALSE)


### 3. GRF ----
cf <- grf::causal_forest(x, y, w)
rf <- grf::regression_forest(x, y)

# prepare a proxy object
grf.obj <- list()
grf.obj$risk.baseline <- as.numeric(rf$predictions)
grf.obj$predicted.absolute.benefit <- as.numeric(cf$predictions)
grf.obj$inputs <- list(X = x, w = w, y = y)
pdf(file = paste0(getwd(), "/tester-plots/const-cf-calibration-absolute.pdf"))
calibration.plot(grf.obj, relative = FALSE, title = "Causal Forest, Calibration: Predicted Absolute Benefit") # TODO: make doubly robust
dev.off()
grf::average_treatment_effect(cf)


