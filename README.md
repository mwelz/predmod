# predictive-modeling
Predictive modeling of heterogeneous treatment effects in medicine.

# Todo list
- [x] Replace penalized regression in 2nd stage risk model with non-penalized one (no need for penalty in 2nd stage). Make sure to report estimates for standard errors in output (same in effect model).
- [ ] How to do risk prediction in Cox PH models (i.e. predicting $Pr(Y = 1 | X)$)? In old version, it was done as follows: `1 - exp(cumulative_base_hazard)^(exp(lp))`, where `cumulative_base_hazard` is taken from `hdnom::glmnet_basesurv` and `lp` is taken from a separate `glm` Cox PH fit. This is how it is done now in Cox risk modeling and Cox baseline modeling. However, based on the `predict.coxph` documentation from `survival`, in Cox effect modeling, I use `1 - exp(predict(model, type = "expected")`. We need to be consistent. Which one shall we use?
- [ ] Computation of C index outcome in Cox risk, stage 2. In stage 1, we performed maximization of the concordance statistic and the lambda corresponding to the largest cross-validated C statistic was the final lambda. Thus, for stage 2, we need to find the concordance associated with model in stage 2. Check how that works; maybe this helps: https://cran.r-project.org/web/packages/survival/vignettes/concordance.pdf. A similar issue will pop up iin Cox effect modeling.
- [x] If constant relative treatment effect is assumed in risk models, glmnet will throw an error because it requires a matrix with at least two variables. Hence, don't use penalized regression in this case.
- [x] Implement Chernozhukov models
- [x] Replace Poisson models with CoxPH
- [x] Update the analyzer of the predictive models 
- [x] Check calculation of C for benefit; requires matching by predicted benefit
