


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
  
  
  #Baseline estimate Poisson model
  Poissonmodel_Baseline <- glm(y ~ x,offset = log(LY),
                               family = poisson(link = "log"))    
  summary(Poissonmodel_Baseline)
  
  #Check values manually, to double check predictions
  intercept = Poissonmodel_Baseline$coefficients[[1]]
  b1 = Poissonmodel_Baseline$coefficients[[2]]
  b2 = Poissonmodel_Baseline$coefficients[[3]]
  b3 = Poissonmodel_Baseline$coefficients[[4]]
  b4 = Poissonmodel_Baseline$coefficients[[5]]
  b5 = Poissonmodel_Baseline$coefficients[[6]]
  
  x1=x[1,1]
  x2=x[1,2]
  x3=x[1,3]
  x4=x[1,4]
  x5=x[1,5]
  LY1 = LY[1]
  
  #Linear predictor poisson model
  LP_Poisson_person_1= intercept+(b1*x1)+(b2*x2)+(b3*x3)+(b4*x4)+(b5*x5) #-2.586006
  #However, Poisson model assumes log(E(Y|x))=(a+Bx). So E(Y|X)=exp(a+Bx)
  Exp_LP_Poisson_person_1 = exp(LP_Poisson_person_1) #0.07532028 ; prediction per personyear
  Exp_LP_Poisson_person_1_LY = Exp_LP_Poisson_person_1*LY1 #0.2863868 ; prediction given observed number of personyears
  
  #Compare with predict() function
  Predict_values_baseline = predict.glm(Poissonmodel_Baseline, type = "response")
  Predict_values_baseline[[1]] #0.2863868 #predictions take into account personyears
  
  LP_Poisson_estimates = intercept+(b1*x[,1])+(b2*x[,2])+(b3*x[,3])+(b4*x[,4])+(b5*x[,5])
  Exp_LP_Poisson_estimates = exp(intercept+(b1*x[,1])+(b2*x[,2]+(b3*x[,3])+(b4*x[,4])+(b5*x[,5])))
  
  
  
  
  #2nd stage
  Poissonmodel_stage2 <- glm(y ~ w,offset = (LP_Poisson_estimates),
                             family = poisson(link = "log")) 
  summary(Poissonmodel_stage2 )
  
  #Coefficients:
  #  Estimate Std. Error z value Pr(>|z|)    
  #(Intercept)  1.75722    0.01936   90.75   <2e-16 ***
  #  w           -0.35935    0.02991  -12.01   <2e-16 ***
  
  #Check how using exponentiated values affect estimates
  Poissonmodel_stage2_exp <- glm(y ~ w,offset = (Exp_LP_Poisson_estimates),
                                 family = poisson(link = "log")) 
  summary(Poissonmodel_stage2_exp )
  
  #Estimate Std. Error z value Pr(>|z|)    
  #(Intercept) -0.70951    0.01936  -36.64   <2e-16 ***
  #  w           -0.35174    0.02991  -11.76   <2e-16 ***
  
  #mostly afffects value of the intercept
  
  Predictdata_LP= data.frame(x,y,LY,w)
  w.rev           <- ifelse(w == 1, 0, 1)
  Predictdata_LP_wflipped= data.frame(x,y,LY,w=w.rev)
  
  #Predicted values for the treatment Poisson model
  LP_Poisson_w=predict.glm(Poissonmodel_stage2,Predictdata_LP, type = "response")
  LP_Poisson_wflipped=predict.glm(Poissonmodel_stage2,Predictdata_LP_wflipped, type = "response")
  
  test = exp( 1.7572+-0.3594*w[1]+  (LP_Poisson_estimates[1]))
  
  
  LP_Poisson_w[1]
  # absolute predicted benefit
  pred.ben.abs.raw <- LP_Poisson_w - LP_Poisson_wflipped
  pred.ben.abs     <- ifelse(w == 1, -pred.ben.abs.raw, pred.ben.abs.raw)
  
  
  
  # relative predicted benefit
  pred.ben.rel.raw <- LP_Poisson_w / LP_Poisson_wflipped
  pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)
  
  pois.ate.hat = mean(pred.ben.abs) #0.1637023
  rel.hat = mean(pred.ben.rel) #0.6981297
  
  #estimated rates
  estimated_w0_rate_per_1000=sum(LP_Poisson_w[w==0])/sum(LY[w==0])*1000
  estimated_w1_rate_per_1000=sum(LP_Poisson_w[w==1])/sum(LY[w==1])*1000
  
  Estimated_reduction_per_1000 = estimated_w1_rate_per_1000-estimated_w0_rate_per_1000
  
  #Effect modeling
  #Baseline estimate Poisson model
  Poissonmodel_Effect <- glm(y ~ x+w,offset = log(LY),
                             family = poisson(link = "log"))    
  summary(Poissonmodel_Effect)
  
  
  
  Predictdata_LP= data.frame(x,y,LY,w)
  w.rev           <- ifelse(w == 1, 0, 1)
  Predictdata_LP_wflipped= data.frame(x,y,LY,w=w.rev)
  
  #Predicted values for the treatment Poisson model
  Effect_Poisson=predict.glm(Poissonmodel_Effect,Predictdata_LP, type = "response")
  Effect_Poisson_wflipped=predict.glm(Poissonmodel_Effect,Predictdata_LP_wflipped, type = "response")
  
  
  # absolute predicted benefit
  pred.ben.abs.raw_effect <- Effect_Poisson - Effect_Poisson_wflipped
  pred.ben.abs_effect     <- ifelse(w == 1, -pred.ben.abs.raw_effect, pred.ben.abs.raw_effect)
  pois.effect.ate.hat = mean(pred.ben.abs_effect ) #0.1765962
  
  
  # relative predicted benefit
  pred.ben.effect.rel.raw <- Effect_Poisson / Effect_Poisson_wflipped
  pred.ben.effect.rel     <- ifelse(w == 1, pred.ben.effect.rel.raw, 1 / pred.ben.effect.rel.raw)
  rel.effect.hat = mean(pred.ben.effect.rel) #0.6789596
  
  #estimated rates
  estimated_w0_effect_rate_per_1000=sum(Effect_Poisson[w==0])/sum(LY[w==0])*1000
  estimated_w1_effect_rate_per_1000=sum(Effect_Poisson[w==1])/sum(LY[w==1])*1000
  
  Estimated_reduction_effect_per_1000 = estimated_w1_effect_rate_per_1000-estimated_w0_effect_rate_per_1000
  
  
