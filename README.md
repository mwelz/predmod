# predictive-modeling
Predictive modeling of heterogeneous treatment effects in medicine.

# Todo list

- [ ] How to do risk prediction in Cox PH models (i.e. predicting $Pr(Y = 1 | X)$)? In old version, it was done as follows: `1 - exp(cumulative_base_hazard)^(exp(lp))`, where `cumulative_base_hazard` is taken from `hdnom::glmnet_basesurv` and `lp` is taken from a separate `glm` Cox PH fit. This is how it is done now in Cox risk modeling and Cox baseline modeling. However, based on the `predict.coxph` documentation from `survival`, in Cox effect modeling, I use `1 - exp(predict(model, type = "expected")`. We need to be consistent. Which one shall we use?
- [x] Implement Chernozhukov models
- [x] Replace Poisson models with CoxPH
- [ ] Write my own Hmisc::cstat to avoid the stupid package clash
- [ ] Equation 3 in rekkas2019 is discussed in the publication. Implement!
- [x] Update the analyzer of the predictive models 
- [x] Check calculation of C for benefit; requires matching by predicted benefit
