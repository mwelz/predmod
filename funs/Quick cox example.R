
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
# Assume base 5 Lifeyears per person, assuming higher risk for disease is correlated with lower life-years and that those who die from disease lose 10% of their lifeyears left
Base_LY   <- 5
LY_coeffs <- c(0.5, -0.3, 0.7, -0.1, 0.4) #reversed signs of the coeffs for diseaseas placeholder
LY_intermediate <- (Base_LY + x  %*% LY_coeffs) * ifelse(y == 1, 0.9, 1) 
# set LY of 0 or less to 0 LY
LY <- ifelse(LY_intermediate <=0, 0, LY_intermediate) 																																												  
max(LY)

#Take a look at the combined dataset
combined_data=data.frame(x,y,LY,w)

Observed_death_rate_1000_PY_treatment = (sum(y[w==1])/sum(LY[w==1]))*1000 #78.9178 deaths per 1,000 person-years
Observed_death_rate_1000_PY_no_treatment = (sum(y[w==0])/sum(LY[w==0]))*1000 #115.1633 deaths per 1,000 person-years

Observed_relative_reduction = Observed_death_rate_1000_PY_treatment / Observed_death_rate_1000_PY_no_treatment #0.6852685
Observed_absolute_reduction_1000_PY = Observed_death_rate_1000_PY_treatment - Observed_death_rate_1000_PY_no_treatment #-36.24553 or, a -0.03624553 reduction per person-year

predictiontimeframe = 4
alpha=1
lifeyears=LY

X=x
offset.lp = TRUE

#NB: ATE and relative benefits are time dependent
two_year = cox.risk.modeling(x, w, y, lifeyears, predictiontimeframe = 2, alpha, offset.lp = TRUE)
four_year = cox.risk.modeling(x, w, y, lifeyears, predictiontimeframe = 4, alpha, offset.lp = TRUE)
six_year = cox.risk.modeling(x, w, y, lifeyears, predictiontimeframe = 6, alpha, offset.lp = TRUE)
eight_year = cox.risk.modeling(x, w, y, lifeyears, predictiontimeframe = 8, alpha, offset.lp = TRUE)

two_year$ate.hat #0
mean(two_year$predicted.relative.benefit) #NAN 

four_year$ate.hat #-0.006247984
mean(four_year$predicted.relative.benefit) #1.096231

six_year$ate.hat #-0.1015129
mean(six_year$predicted.relative.benefit) # 0.8846003

eight_year$ate.hat #-0.0008124408
mean(eight_year$predicted.relative.benefit) #0.9991255




two_year_e= cox.effect.modeling(x, w, y, LY, predictiontimeframe=2,alpha = 1)
two_year_e$ate.hat #0
mean(two_year_e$predicted.relative.benefit) #NAN 

four_year_e= cox.effect.modeling(x, w, y, LY, predictiontimeframe=4,alpha = 1)
four_year_e$ate.hat # -0.05487191
mean(four_year_e$predicted.relative.benefit) #0.6025049 

six_year_e= cox.effect.modeling(x, w, y, LY, predictiontimeframe=6,alpha = 1)
six_year_e$ate.hat # -0.05487191
mean(six_year_e$predicted.relative.benefit) #0.6025049 