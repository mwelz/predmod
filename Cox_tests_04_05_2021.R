
#designate logistic function
logistic <- function(x) 1 / (1 + exp(-x))

### 0.1. Data generation ---- 10,000 persons, 5 covariates
set.seed(2)
n <- 10000
p <- 5

# treatment assignment, 50% treatment, 50% control.
w <- rbinom(n, 1, 0.5) 

# covariates for outcome variable (y)
x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))

# coefficients (including an intercept)
theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)

# compute Pr(Y = 1 | X) for each individual (with noise)
eps <-  rnorm(n, mean = 0, sd = 0.5)
pi0 <- logistic(as.numeric(cbind(1,x) %*% theta) + eps)

#assume a true relative constant risk reduction of 30%
scaling <- 0.7
pi1 <- pi0 * scaling 

# create binary outcomes
y0 <- rbinom(n, 1, pi0)
y1 <- rbinom(n, 1, pi1)
y  <- ifelse(w == 1, y1, y0) # observed outcome

### 0.2. Generate lifeyears for each individual, assume base 5 lifeyears per person
Base_LY   <- 5
#Generate set of covariates unrelated to original x for lifeyears
x2 <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
#assume that those who die from the disease lose 10% of their lifeyears left
LY_coeffs <- c(-0.3, 0.2, -0.4, 0.2, -0.4) 
LY_intermediate <- (Base_LY + x2  %*% LY_coeffs) * ifelse(y == 1, 0.9, 1) 

# set LY of 0 or less to 0 LY
LY <- ifelse(LY_intermediate <=0, 0, LY_intermediate)
max_lifeyears = max(LY) #7.422949
minimum_lifeyears = min(LY) #2.211346



#prediction timeframe
predictiontimeframe = 4 #assume 4 year timeframe for predictions; assume individuals are censored after 4 years
LY_predictiontimeframe <- ifelse(LY <=predictiontimeframe, LY, predictiontimeframe)
y_predictiontimeframe <- ifelse(LY <=predictiontimeframe, y, 0)

#Observed with general LY
Observed_death_rate_1000_PY_timeframe_treatment = (sum(y_predictiontimeframe[w==1])/sum(LY_predictiontimeframe[w==1]))*1000 #20.175 deaths per 1,000 person-years
Observed_death_rate_1000_PY_timeframe_no_treatment = (sum(y_predictiontimeframe[w==0])/sum(LY_predictiontimeframe[w==0]))*1000 #29.6742 deaths per 1,000 person-years
Observed_absolute_reduction_timeframe_1000_PY = Observed_death_rate_1000_PY_timeframe_treatment - Observed_death_rate_1000_PY_timeframe_no_treatment #-9.499197  per 1,000 person-years

Observed_relative_reduction_timeframe = Observed_death_rate_1000_PY_timeframe_treatment / Observed_death_rate_1000_PY_timeframe_no_treatment #0.6798836



#Estimate models with RMS
wx <- unname(cbind(x, w = ifelse(w == 1, 0, 1)))
glmsurvival.coxeffect.obj <-survival::Surv(LY_predictiontimeframe, y_predictiontimeframe)
colnames(glmsurvival.coxeffect.obj) <- c("time", "status")
final.model.fit <- rms::cph(glmsurvival.coxeffect.obj~ wx,y=TRUE,x=TRUE)
ph_check_rms=survival::cox.zph(final.model.fit , transform="km",terms=FALSE)
print(ph_check_rms) #PH assumption p-values check out
plot(ph_check_rms)

Dxy = final.model.fit$stats[[9]]
abs.Dxy = abs(Dxy)
c.stat=(abs.Dxy/2)+0.5 #0.6254992

#Get linear predictor and probabilities
lp=final.model.fit$linear.predictors
lp_alternate = final.model.fit$coefficients * wx 
difference_lp = lp-lp_alternate #difference in linear predictors; potentially due to normalization in rms function?

#check with general survival model:
res.cox <- survival::coxph(glmsurvival.coxeffect.obj ~ wx)
summary(res.cox)
ph_check_surv=survival::cox.zph(res.cox , transform="km",terms=FALSE)
print(ph_check_surv) #PH assumption p-values check out
plot(ph_check_surv)

lp_surv = res.cox$linear.predictors
difference_lp2 = lp-lp_surv #uniform difference of about ~ -0.2032 for RMS vs general survival package; difference due to centering constant?

#expected from res.cox
surv_lp= predict(res.cox, type="lp")
surv_exp= predict(res.cox, type="expected")

#Check E/O calibration from res.cox
sum(surv_exp[w==1]) #404
sum(y_predictiontimeframe[w==1]) #404
sum(surv_exp[w==0]) #579
sum(y_predictiontimeframe[w==0]) #579
sum(surv_exp) #983
sum(y_predictiontimeframe)#983



#Manual calculations of absolute risk: potentially off due to normalization and or centering?
#Manual calculation of absolute risk; rms
basehazard_rms <-  survival::basehaz(final.model.fit,centered=FALSE) #assume centered = false, but no real difference if set to TRUE
Hazard_t_index_rms = match(TRUE, round(basehazard_rms$time,2) == predictiontimeframe, nomatch = NA)
probs_rms  =  1-exp(-basehazard_rms[Hazard_t_index_rms,1])^exp(lp) 
mean(abs(probs_rms-surv_exp)) #0.02633046 ; wrong, likely in part due to lp
max(abs(probs_rms-surv_exp)) #0.1836317

#Manual calculation of absolute risk; general survival
basehazard_surv <-  survival::basehaz(res.cox,centered=FALSE)
Hazard_t_index_surv  = match(TRUE, round(basehazard_surv$time,2) == predictiontimeframe, nomatch = NA)
probs_surv   =  1-exp(-basehazard_surv[Hazard_t_index_surv,1])^exp(surv_lp) 
mean(abs(probs_surv-surv_exp)) #0.01345942 ; wrong, likely in part due to lp
max(abs(probs_surv-surv_exp)) #0.2183548



#glmnet
glmsurvival.coxeffect.obj <-survival::Surv(LY_predictiontimeframe, y_predictiontimeframe)
colnames(glmsurvival.coxeffect.obj) <- c("time", "status")
glm.cox <- glmnet::cv.glmnet(wx, glmsurvival.coxeffect.obj, family = "cox", type.measure = "C", alpha = 1)
coefs.glm.cox     <- glmnet::coef.glmnet(glm.cox, s = "lambda.min")
lambda        <- glm.cox$lambda.min
lp.glm <-predict(glm.cox$glmnet.fit, s=lambda, wx,type="link") #linear predictors?

