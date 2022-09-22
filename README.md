# predictive-modeling
Predictive modeling of heterogeneous treatment effects in medicine.

# Todo list
- [ ] Can we bootstrap standard error of relative risk estimate? No we cannot, so don't compute SE for relative risks
- [ ] WA audio
- [ ] external validation function: predict methods!

# Examples
## Complete case

```R
library(predmod)

set.seed(1)
lung <- survival::lung
lung <- na.omit(lung)
time <- lung$time
status <- lung$status - 1L
X <- as.matrix(lung[, 4:10]) # randomly assign treatment status
W <- rbinom(nrow(X), 1, 0.5)

## 1. cross-sectional models (disregard time)
# fit
crss_risk <- risk_model(X = X, status = status, w = W)
crss_effect <- effect_model(X = X, status = status, w = W, interacted = c("sex", "age"))
crss_grf <- grf_model(X = X, status = status, w = W, num_trees = 500)

# plot
calibration_plot(crss_risk)
calibration_plot(crss_effect)
calibration_plot_grf(crss_grf)

## 2. survival models
# fit
surv_risk   <- risk_model_survival(X = X, status = status, time = time, w = W)
surv_effect <- effect_model_survival(X = X, status = status, time = time, w = W, interacted = c("age", "sex"))
surv_grf    <- grf_model_survival(X = X, status = status, time = time, w = W, num_trees = 500)

# plot
calibration_plot(surv_risk)
calibration_plot(surv_effect)
calibration_plot_grf(surv_grf, baseline_risk = surv_risk$risk$baseline)

```

## Imputation 
