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
#' Authors: kth, mwelz
#' Last changed: Feb 10, 2021, by mwelz   
################################################################################
rm(list = ls()) ; cat("\014")

# load the helper functions
source(paste0(getwd(), "/funs/estimation-funs.R"))

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
risk.model$c.index # 0.713

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
effect.model$c.index # 0.712

### 3. GRF ----
grf.obj <- grf.modeling(X = x, y = y, w = w)
grf.obj$ate # -0.166
ate.ci.lo <- grf.obj$ate.hat - grf.obj$ate.hat.se * qt(0.975, df = n - p)
ate.ci.up <- grf.obj$ate.hat + grf.obj$ate.hat.se * qt(0.975, df = n - p)
(ate.ci.lo <= ate) & (ate <= ate.ci.up) # ATE is in 95% CI

grf.obj$c.index # 0.679 # but this is just experimental!
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


### 5.0 Poisson regression ----
#load rateratio.test
library(rms)

#Regular glm


#Glm (note capital G) from rms; need to check what causes differences
#Poissonmodel2 <- Glm(y ~ x+w,
#                    offset = log(LY),
#                    family = poisson(link = "log"))    
#Poissonmodel2

Poissonmodel_Baseline <- glm(y ~ x,offset = log(LY),
                          family = poisson(link = "log"))    
summary(Poissonmodel_Baseline )

#Values for x coefficients are off; in part affected by LY offset (see further analyses below)



#Make a separate dataframe to allow for predictions from the Poisson model
Predictdata_Baseline = data.frame(x,y,LY)

#Predicted values for the Poisson model
LP_Poisson=predict(Poissonmodel_Baseline,name=Predictdata_Baseline, type = "response")

Predictdata_LP= data.frame(x,y,LY,LP_Poisson,w)
w.rev           <- ifelse(w == 1, 0, 1)
Predictdata_LP_wflipped= data.frame(x,y,LY,LP_Poisson,w=w.rev)

#Not yet sure how to do the offset for lifeyears here; the LP does take those into account though.
Poissonmodel_Treat <- glm(y ~ w,offset = log(LP_Poisson),
                   family = poisson(link = "log")) 
summary(Poissonmodel_Treat )



#Predicted values for the treatment Poisson model
LP_Poisson=predict.glm(Poissonmodel_Treat,Predictdata_LP, type = "response")
LP_Poisson_wflipped=predict.glm(Poissonmodel_Treat,Predictdata_LP_wflipped, type = "response")


# absolute predicted benefit
pred.ben.abs.raw <- LP_Poisson - LP_Poisson_wflipped
pred.ben.abs     <- ifelse(w == 1, -pred.ben.abs.raw, pred.ben.abs.raw)

# relative predicted benefit
pred.ben.rel.raw <- LP_Poisson / LP_Poisson_wflipped
pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)

pois.ate.hat = mean(pred.ben.abs)
#

Predictdata_LP= data.frame(x,y,LY,LP_Poisson,w)
Treatmentsubset=subset(Predictdata_LP,w==1)
No_Treatmentsubset=subset(Predictdata_LP,w==0)

#rates per 1,000 personyears
Rate_treatment_1000 =  (sum(Treatmentsubset$y)/sum(Treatmentsubset$LY))*1000
Rate_no_treatment_1000 = (sum(No_Treatmentsubset$y)/sum(No_Treatmentsubset$LY))*1000

rel.hat = mean(pred.ben.rel)
rel.hat2 = Rate_treatment_1000/Rate_no_treatment_1000

##No LY Offset analyses

Poissonmodel_Baseline_no_offset <- glm(y ~ x,
                             family = poisson(link = "log"))    
summary(Poissonmodel_Baseline_no_offset )









### 6.0 Cox proportional hazard----
S <- survival::Surv(time = LY, event = y) #Survival time until event; if 0, right censored.
f <- rms::cph(S ~ x + w, surv = TRUE)
g <- rms::Survival(f)   # g is a function

ftwo <- survival::coxph(S ~ x + w)

predict(ftwo,type="expected")

# TODO: c-statistic
# TODO: make calibration plot doubly robust when we use causal forest




