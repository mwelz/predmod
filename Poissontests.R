
Poisson_runs <- function(startinglifeyears,diseasereduction,relative_effect) 
{
  
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
scaling <- relative_effect
pi1 <- pi0 * scaling 
hte <- pi1 - pi0

# true ATE
ate <- mean(pi1) - mean(pi0)

# create binary outcomes
y0 <- rbinom(n, 1, pi0)
y1 <- rbinom(n, 1, pi1)
y  <- ifelse(w == 1, y1, y0) # observed outcome



### 0.2. Generating lifeyears
# Assume base 5 Lifeyears per person, assuming higher risk for disease is correlated with lower life-years and that those who die from disease lose 10% of their lifeyears left
Base_LY   <- startinglifeyears
LY_coeffs <- c(0.5, -0.3, 0.7, -0.1, 0.4) #reversed signs of the coeffs for diseaseas placeholder
LY_intermediate <- (Base_LY + x  %*% LY_coeffs) * ifelse(y == 1, diseasereduction, 1) 
# set LY of 0 or less to 0 LY
LY <- ifelse(LY_intermediate <=0, 0, LY_intermediate) 																																												  



#Baseline estimate Poisson model
Poissonmodel_Baseline <- glm(y ~ x,offset = log(LY),
                             family = poisson(link = "log"))    
 
#Check values manually, to double check predictions
intercept = Poissonmodel_Baseline$coefficients[[1]]
b1 = Poissonmodel_Baseline$coefficients[[2]]
b2 = Poissonmodel_Baseline$coefficients[[3]]
b3 = Poissonmodel_Baseline$coefficients[[4]]
b4 = Poissonmodel_Baseline$coefficients[[5]]
b5 = Poissonmodel_Baseline$coefficients[[6]]


LP_Poisson <-intercept+(b1*x[,1])+(b2*x[,2])+(b3*x[,3])+(b4*x[,4])+(b5*x[,5])



#2nd stage
Poissonmodel_stage2 <- glm(y ~ w,offset = LP_Poisson, family = poisson(link = "log")) 



Predictdata_LP = data.frame(x,y,LY,w)
w.rev           <- ifelse(w == 1, 0, 1)
Predictdata_LP_wflipped = data.frame(x,y,LY,w=w.rev)

#Predicted values for the treatment Poisson model
#Poisson_w = predict.glm(Poissonmodel_stage2,Predictdata_LP, type = "response")
#Poisson_wflipped = predict.glm(Poissonmodel_stage2,Predictdata_LP_wflipped, type = "response")

Poisson_w = exp(Poissonmodel_stage2$coefficients[[1]]+Poissonmodel_stage2$coefficients[[2]]*w + (LP_Poisson))
Poisson_wflipped = exp(Poissonmodel_stage2$coefficients[[1]]+Poissonmodel_stage2$coefficients[[2]]*w.rev +  (LP_Poisson))

# absolute predicted benefit
pred.ben.abs.raw <- Poisson_w - Poisson_wflipped
pred.ben.abs     <- ifelse(w == 1, -pred.ben.abs.raw, pred.ben.abs.raw)



# relative predicted benefit
pred.ben.rel.raw <- Poisson_w / Poisson_wflipped
pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)

pois.ate.hat = mean(pred.ben.abs) 
rel.hat = mean(pred.ben.rel) 

#estimated rates
estimated_w0_rate_per_1000=sum(Poisson_w[w==0])/sum(LY[w==0])*1000
estimated_w1_rate_per_1000=sum(Poisson_w[w==1])/sum(LY[w==1])*1000

Estimated_reduction_per_1000 = estimated_w1_rate_per_1000-estimated_w0_rate_per_1000

#Effect modeling
#Baseline estimate Poisson model
Poissonmodel_Effect <- glm(y ~ x+w,offset = log(LY),
                             family = poisson(link = "log"))    




Predictdata_LP= data.frame(x,y,LY,w)
w.rev           <- ifelse(w == 1, 0, 1)
Predictdata_LP_wflipped= data.frame(x,y,LY,w=w.rev)

#Predicted values for the treatment Poisson model
Effect_Poisson=predict.glm(Poissonmodel_Effect,Predictdata_LP, type = "response")
Effect_Poisson_wflipped=predict.glm(Poissonmodel_Effect,Predictdata_LP_wflipped, type = "response")


# absolute predicted benefit
pred.ben.abs.raw_effect <- Effect_Poisson - Effect_Poisson_wflipped
pred.ben.abs_effect     <- ifelse(w == 1, -pred.ben.abs.raw_effect, pred.ben.abs.raw_effect)
pois.effect.ate.hat = mean(pred.ben.abs_effect ) 


# relative predicted benefit
pred.ben.effect.rel.raw <- Effect_Poisson / Effect_Poisson_wflipped
pred.ben.effect.rel     <- ifelse(w == 1, pred.ben.effect.rel.raw, 1 / pred.ben.effect.rel.raw)
rel.effect.hat = mean(pred.ben.effect.rel) 
#estimated rates
estimated_w0_effect_rate_per_1000=sum(Effect_Poisson[w==0])/sum(LY[w==0])*1000
estimated_w1_effect_rate_per_1000=sum(Effect_Poisson[w==1])/sum(LY[w==1])*1000

Estimated_reduction_effect_per_1000 = estimated_w1_effect_rate_per_1000-estimated_w0_effect_rate_per_1000


observed_rate_per_1000_treat=sum(y[w==0])/sum(LY[w==0])*1000
observed_rate_per_1000_no_treat=sum(y[w==1])/sum(LY[w==1])*1000
Observed_reduction_per_1000 = observed_rate_per_1000_no_treat-observed_rate_per_1000_treat

return(list(
  generated.ate = ate,
  generated.rel = relative_effect,
  risk.ate = pois.ate.hat,
  risk.rel = rel.hat,
  effect.ate = pois.effect.ate.hat,
  effect.rel = rel.effect.hat,
  generated.reduction.rate.1000 = Observed_reduction_per_1000,
  estimated.reduction.risk.rate.1000 =Estimated_reduction_per_1000,
  estimated.reduction.effect.rate.1000 = Estimated_reduction_effect_per_1000)
)

}

