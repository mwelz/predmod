rm(list = ls()) ; cat("\014")

# get data
library(survival)
lung <- survival::lung

# 1 means having the event
lung$status <- ifelse(lung$status == 1, 1, 0)

# for now, kick out duplicate times
lung <- lung[!duplicated(lung$time),]

# run model
res.cox <- coxph(Surv(time, status) ~ sex + age, data = lung )
res.cox

x = cbind(sex = lung$sex, age = lung$age)
x = as.matrix(x) # x needs to be a matrix!
time = lung$time
beta = c(-0.1, 0.1)
status = lung$status
i=1
Ri <- which(time >= time[i])

# x needs to be a matrix with more than one column, otherwise error due

sum_over_Ri <- function(Ri, x, beta){
  
  sum(sapply(Ri, function(j) as.numeric(exp(t(x[j,]) %*% beta))))
  
}


sum_over_Ri_gradient <- function(Ri, x, beta){
  
  rowSums(sapply(Ri, function(j) x[j,] * as.numeric(exp(t(x[j,]) %*% beta))))
  
}


logL_individual_contribution <- function(x, i, beta, time){

  Ri <- which(time >= time[i])
  
  as.numeric(t(x[i,]) %*% beta - log(sum_over_Ri(Ri, x, beta)))
  
} # FUN


logL_individual_contribution_gradient <- function(x, i, beta, time){
  
  Ri <- which(time >= time[i])
  
  as.numeric(x[i,] - sum_over_Ri_gradient(Ri, x, beta) / sum_over_Ri(Ri, x, beta))
  
} # FUN


neglogL_gradient <- function(x, beta, time, status){
  
  failures <- which(status == 1)
  -rowSums(
    sapply(failures, function(i)  logL_individual_contribution_gradient(x, i, beta, time))
    )

} 


neglogL <- function(x, beta, time, status){
  
  failures <- which(status == 1)
  -sum(sapply(failures, function(i)  logL_individual_contribution(x, i, beta, time) ))
  
}

#neglogL(x, beta, time, status)
#neglogL_gradient(x, beta, time, status)

foo = optim(par = rep(0, ncol(x)),
            fn = neglogL, gr = neglogL_gradient, method = "BFGS",
            x = x, time = time, status = status)
foo$par
res.cox$coefficients
