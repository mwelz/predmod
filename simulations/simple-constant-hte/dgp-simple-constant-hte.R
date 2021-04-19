################################################################################
#' dgp-simple-constant-hte
#' 
#' In this script, we employ a simulation design to model the most simple DGP.
#' We test the predictive models in this simple situation.
#' The DGP has the following properties:
#' 
#' - it is an RCT
#' - 5 iid normally distributed covariates
#' - constant relative risk reduction of 30%. That is, 
#'   Pr(Y = 1 | X, W = 1) = 0.7 * Pr(Y = 1 | X, W = 0),
#'   where Y = 1 indicates a death.
#' - lifeyears are also generated (needed for duration models)
#' 
#' Authors: mwelz, kth
#' Last changed: Feb 17, 2021, by mwelz   
################################################################################
rm(list = ls()) ; cat("\014")

# load the helper functions
source(paste0(getwd(), "/funs/estimation-funs.R"))
source(paste0(getwd(), "/funs/Poiss_test_estimation-funs.R")) 
source(paste0(getwd(), "//funs/Poiss_treat_test_estimation-funs.R")) 
source(paste0(getwd(), "/funs/cox_risk_modelling.R")) 
source(paste0(getwd(), "//funs/cox_effect_modelling.R")) 
# define the logistic function
logistic <- function(x) 1 / (1 + exp(-x))

### 0.1. Data generation ----
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


### 0.2. Generating lifeyears
# Assume base 5 Lifeyears per person, assuming higher risk for lung cancer is correlated with lower life-years and that those who die from lung cancer lose 10% of their lifeyears left
Base_LY   <- 5
LY_coeffs <- c(0.5, -0.3, 0.7, -0.1, 0.4) #reversed signs of the coeffs for lung cancer as placeholder
LY_intermediate <- (Base_LY + x  %*% LY_coeffs) * ifelse(y == 1, 0.9, 1) 
# set LY of 0 or less to 0 LY
LY <- ifelse(LY_intermediate <=0, 0, LY_intermediate) 																																												  

### 1. risk modeling ----
risk.model <- risk.modeling(X = x, w = w, y = y, alpha = 1, offset.lp = TRUE)

# # make plots (currently commented out) 
# pdf(file = paste0(getwd(), "/plots/const-rm-calibration-relative.pdf"))
# calibration.plot(risk.model, relative = TRUE, title = "Risk Model, Calibration: Predicted Relative Benefit")
# dev.off()
# 
# pdf(file = paste0(getwd(), "/plots/const-rm-calibration-absolute.pdf"))
# calibration.plot(risk.model, relative = FALSE, title = "Risk Model, Calibration: Predicted Absolute Benefit")
# dev.off()
# 
# pdf(file = paste0(getwd(), "/plots/const-rm-subgroup-absolute.pdf"))
# subgroup.plot(risk.model, x[,1], relative = FALSE)
# dev.off()
# 
# pdf(file = paste0(getwd(), "/plots/const-rm-subgroup-relative.pdf"))
# subgroup.plot(risk.model, x[,1], relative = TRUE)
# dev.off()

risk.model$ate.hat # 0.164 (average of predicted abolute benefits)
mean(risk.model$predicted.relative.benefit) # 0.68
risk.model$c.index.benefit #0.4555163

### 2. effect modeling ----
effect.model <- effect.modeling(X = x, w = w, y = y, alpha = 1) 
# 
# # make plots (currently commented out)
# calibration.plot(effect.model, relative = TRUE, title = "Effect Model, Calibration: Predicted Relative Benefit")
# calibration.plot(effect.model, relative = FALSE, title = "Effect Model, Calibration: Predicted Absolute Benefit")
# 
# subgroup.plot(effect.model, x[,1], relative = TRUE)
# subgroup.plot(effect.model, x[,1], relative = FALSE)

effect.model$ate.hat # 0.161 (average of predicted abolute benefits)
mean(effect.model$predicted.relative.benefit) # 0.675
effect.model$c.index.benefit #0.4676652

### 3. GRF ----
grf.obj <- grf.modeling(X = x, y = y, w = w)
grf.obj$ate.hat # -0.166
ate.ci.lo <- grf.obj$ate.hat - grf.obj$ate.hat.se * qt(0.975, df = n - p)
ate.ci.up <- grf.obj$ate.hat + grf.obj$ate.hat.se * qt(0.975, df = n - p)
(ate.ci.lo <= ate) & (ate <= ate.ci.up) # ATE is in 95% CI

grf.obj$c.index.benefit # 0.4798283
grf.calibration <- calibration.plot.grf(grf.obj)


# cf <- grf::causal_forest(x, y, w)
# rf <- grf::regression_forest(x, y)
# 
# # prepare a proxy object as input for the calibration plots
# grf.obj <- list()
# grf.obj$risk.baseline <- as.numeric(rf$predictions)
# grf.obj$predicted.absolute.benefit <- as.numeric(cf$predictions)
# grf.obj$inputs <- list(X = x, w = w, y = y)
# grf::average_treatment_effect(cf) # -0.165. Disadvantage: No estimation of relative risk feasible.

# # calibration plot (commented out)
# pdf(file = paste0(getwd(), "/plots/const-cf-calibration-absolute.pdf"))
# calibration.plot(grf.obj, relative = FALSE, title = "Causal Forest, Calibration: Predicted Absolute Benefit") # TODO: make doubly robust
# dev.off()

### 4.0 Rate-ratio ----

## overall Rate-Ratio
overall.rateratio <- rate.ratio(y = y, w = w, lifeyears = LY)$rate.ratio # 0.6852685; close to 0.7. But, will  be more/less favourable depending on how lifeyears are affected

## Example: those who die of lung cancer lose 25% of their remaining lifeyears instead of 0.9
LY_intermediate2 <- (Base_LY + x  %*% LY_coeffs) * ifelse(y == 1, 0.75, 1) 
LY2 <- ifelse(LY_intermediate2 <=0, 0, LY_intermediate2) # set LY of 0 or less to 0 LY
overall.rateratio2 <- rate.ratio(y = y, w = w, lifeyears = LY2)$rate.ratio # 0.6646533 


## Examples for rate-ratios for subgroups:
# Rate-ratio for an example subgroup: x1 < 0 
rate.ratio(y = y, w = w, lifeyears = LY, subgroup = x[,1] < 0)$rate.ratio # 0.6771719

# Rate-ratio for an example subgroup: x1 >= 0 
rate.ratio(y = y, w = w, lifeyears = LY, subgroup = x[,1] >= 0)$rate.ratio # 0.6862652


### 5.0 Cox regression ----

#NB: benefits are time dependent
plot(survival::survfit(survival::Surv(LY, y) ~ w))

#risk
two_year_cox_risk = cox.risk.modeling(X=x, w, y, lifeyears=LY, predictiontimeframe = 2, alpha, offset.lp = TRUE)
two_year_cox_risk$ate.hat #0
two_year_cox_risk$c.index.benefit #0.5
mean(two_year_cox_risk$predicted.relative.benefit)#NaN

four_year_cox_risk = cox.risk.modeling(X=x, w, y, lifeyears=LY, predictiontimeframe = 4, alpha, offset.lp = TRUE)
four_year_cox_risk$ate.hat #-0.006247984
four_year_cox_risk$c.index.benefit #0.4372186
mean(four_year_cox_risk$predicted.relative.benefit)#1.099922

six_year_cox_risk = cox.risk.modeling(X=x, w, y, lifeyears=LY, predictiontimeframe = 6, alpha, offset.lp = TRUE)
six_year_cox_risk$ate.hat #-0.1013502
six_year_cox_risk$c.index.benefit #0.5183192
mean(six_year_cox_risk$predicted.relative.benefit)#0.8848779


eight_year_cox_risk = cox.risk.modeling(X=x, w, y, lifeyears=LY, predictiontimeframe = 8, alpha, offset.lp = TRUE)
eight_year_cox_risk$ate.hat #-0.0008124408
eight_year_cox_risk$c.index.benefit # 0.4974462
mean(eight_year_cox_risk$predicted.relative.benefit)#0.9991255

#cox effect
two_year_cox_effect = cox.effect.modeling(X=x, w, y, lifeyears=LY, predictiontimeframe=2,alpha = 1)
two_year_cox_effect$ate.hat #0
two_year_cox_effect$c.index.benefit #0.5
mean(two_year_cox_effect$predicted.relative.benefit)#NaN

four_year_cox_effect =  cox.effect.modeling(X=x, w, y, lifeyears=LY, predictiontimeframe=4,alpha = 1)
four_year_cox_effect$ate.hat # -0.0548291
four_year_cox_effect$c.index.benefit #0.5371681
mean(four_year_cox_effect$predicted.relative.benefit)#0.6023313

six_year_cox_effect =  cox.effect.modeling(X=x, w, y, lifeyears=LY, predictiontimeframe=6,alpha = 1)
six_year_cox_effect$ate.hat #-0.05720218
six_year_cox_effect$c.index.benefit #0.4582583
mean(six_year_cox_effect$predicted.relative.benefit)#0.9120426


eight_year_cox_effect =  cox.effect.modeling(X=x, w, y, lifeyears=LY, predictiontimeframe=8,alpha = 1)
eight_year_cox_effect$ate.hat #-0.001114145
eight_year_cox_effect$c.index.benefit # 0.4751412
mean(eight_year_cox_effect$predicted.relative.benefit)#0.9987149
