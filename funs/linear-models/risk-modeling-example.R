rm(list = ls()) ; cat("\014")

# load functions
source(paste0(getwd(), "/funs/linear-models/risk-modeling.R"))
source(paste0(getwd(), "/funs/plotmakers/plotmakers.R"))


# make data
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

## test the model
stage1 <- baseline.risk(X, y)
stage2 <- risk.model.stage2(linear.predictor = stage1$linear.predictor, y = y, w = w,
                            z = "linear.predictor", 
                            constant.treatment.effect = FALSE, 
                            intercept = TRUE)

RM <- risk.modeling(X, y, w, alpha = 1)
RM$models$coefficients.stage2

summary(stage2$mod.stage2)$coefficients
stage1$coefficients

# stage2$risk.regular.w
# stage2$risk.flipped.w
# 
# # test the risk model
# alpha <- 1
# intercept.stage.2 <- intercept <- FALSE
# z <- "linear.predictor"
# prediction.timeframe = NULL
# lifeyears = NULL
# constant.treatment.effect=F
# stage1 <- baseline.risk(X, y)
# linear.predictor = stage1$linear.predictor
# 
# RM <- risk.modeling(X = X, y = y, w = w, alpha = alpha, 
#                     intercept.stage.2 = intercept.stage.2, z = z, 
#                     lifeyears = lifeyears, prediction.timeframe = prediction.timeframe)
# 
# calibration.plot(RM)
# subgroup.plot(RM, X[,1])
# 
# RM$models$coefficients.stage2
# RM$average.treatment.effect
# RM$C.statistics
