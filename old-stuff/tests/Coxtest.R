
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


#Take a look at the combined dataset
combined_data=data.frame(x,y,LY,w)

Observed_death_rate_1000_PY_treatment = (sum(y[w==1])/sum(LY[w==1]))*1000 #78.9178 deaths per 1,000 person-years
Observed_death_rate_1000_PY_no_treatment = (sum(y[w==0])/sum(LY[w==0]))*1000 #115.1633 deaths per 1,000 person-years

Observed_relative_reduction = Observed_death_rate_1000_PY_treatment / Observed_death_rate_1000_PY_no_treatment #0.6852685
Observed_absolute_reduction_1000_PY = Observed_death_rate_1000_PY_treatment - Observed_death_rate_1000_PY_no_treatment #-36.24553 or, a -0.03624553 reduction per person-year

Average_PY = mean(LY)                     #4.753821
Average_PY_treatment = mean(LY[w==1])     #4.81719
Average_PY_no_treatment = mean(LY[w==0])  #4.688887



#observed reduction per person-year * average number of personyears in the treatment arm
Average_PY_treatment*(Observed_absolute_reduction_1000_PY/1000) #-0.1746016

survivalobject <- survival::Surv(LY, y)
res.cox <- rms::cph(survivalobject~ x)

x1=x[,1]
x2=x[,2]
x3=x[,3]
x4=x[,4]
x5=x[,5]

res.cox2 <- rms::cph(survivalobject~ x1+x2+x3+x4+x5)
#Get Dxy (Somers' D); Dxy = 2*(c- 0.5) 
#Dxy can be negative in some cases: https://stat.ethz.ch/pipermail/r-help/2011-February/269589.html
#"When predicting relative log hazard, high hazard means short survival time so Dxy is negative."
Dxy = res.cox$stats[[9]]
abs.Dxy = abs(Dxy)
c.stat=(abs.Dxy/2)+0.5

#Proportional hazard assumption testing: need to unpack individual x variables
cox.zph(res.cox , transform="km")
cox.zph(res.cox2 , transform="km")

# glmnet seems to require explicit columnnames for a surival object?
glmsurvival.obj <-Surv(LY, y)
colnames(glmsurvival.obj) <- c("time", "status")
cvfit <- glmnet::cv.glmnet(x, glmsurvival.obj, family = "cox", type.measure = "C", alpha = 1)
coef(cvfit, s = "lambda.min")
#C-stat belonging to lambda.min
lambda.min.index = which(cvfit[["lambda"]] %in% cvfit$lambda.min)
c.stat.lambda.min = cvfit[["cvm"]][lambda.min.index]
