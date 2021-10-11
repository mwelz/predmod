rm(list = ls()) ; cat("\014")
library(cmprsk)

# simulated data to test
set.seed(10)
time <- rexp(200)
status <- sample(0:2,200,replace=TRUE)
x <- matrix(runif(600),nrow=200)
dimnames(x)[[2]] <- c('x1','x2','x3')

#print(z <- crr(ftime,fstatus,cov))
#summary(z)
#z.p <- predict(z,rbind(c(.1,.5,.8),c(.1,.5,.2)))
#plot(z.p,lty=1,color=2:3)

################################################

# Status is binary
km_fit <- function(time, status){
  
  survival::survfit(survival::Surv(time, status) ~ 1, 
                    data = data.frame(time, status),
                    type = "kaplan-meier")
} # FUN

km_predict_survival <- function(x, time){
  
  # x is returned by km_fit()
  summary(x, time = time)$surv
  
} # FUN


# i is in I_k, j is in R_{i,k}. 
Zw <- function(time_i, time_j, km_fit_object, status_j, k){
  
  if((time_j < time_i) & (0 < status_j) & (status_j != k)){
    km_predict_survival(km_fit_object, time_i) / 
      km_predict_survival(km_fit_object, time_j)
  } else 1
  
} # FUN


sum_over_Rik <- function(Rik, i, x, beta, time, status, km_fit_object, k){
  
  sum(sapply(Rik, function(j){ 
    
    as.numeric( Zw(time_i = time[i], time_j = time[j], 
                   km_fit_object = km_fit_object, 
                   status_j = status[j], k = k) * exp(t(x[j,]) %*% beta))
    
    }))
} # FUN


sum_over_Rik_gradient <- function(Rik, i, x, beta, time, status, km_fit_object, k){
  
  rowSums(sapply(Rik, function(j){ 
    
    x[j,] * 
    as.numeric( Zw(time_i = time[i], time_j = time[j], 
                   km_fit_object = km_fit_object, 
                   status_j = status[j], k = k) * exp(t(x[j,]) %*% beta))
    
  }))
} # FUN


logL_individual_contribution_k <- function(i, x, beta, time, status, km_fit_object, k){
  
  Rik <- which(time >= time[i] | (0 < status & status != k))
  
  as.numeric(t(x[i,]) %*% beta - 
               log(sum_over_Rik(Rik, i, x, beta, time, status, km_fit_object, k)))
  
} # FUN



logL_individual_contribution_k_gradient <- function(i, x, beta, time, status, km_fit_object, k){
  
  Rik <- which(time >= time[i] | (0 < status & status != k))
  
  as.numeric(x[i,] - 
               sum_over_Rik_gradient(Rik, i, x, beta, time, status, km_fit_object, k) /
               sum_over_Rik(Rik, i, x, beta, time, status, km_fit_object, k))
  
} # FUN


neglogL_k <- function(x, beta, time, status, k){
  
  delta <- 1 * (status > 0)
  km_fit_object <- km_fit(time = time, status = 1 - delta)
  
  failures.k <- which(status == k)
  -sum(sapply(failures.k, function(i){
    logL_individual_contribution_k(i, x, beta, time, status, km_fit_object, k)
    }))
  
} # FUN

neglogL_k_gradient <- function(x, beta, time, status, k){
  
  delta <- 1 * (status > 0)
  km_fit_object <- km_fit(time = time, status = 1 - delta)
  
  failures.k <- which(status == k)
  
  -rowSums(
    sapply(failures.k, function(i){
      logL_individual_contribution_k_gradient(i, x, beta, time, status, km_fit_object, k)})
  )
  
} # FUN


foo = optim(par = rep(0, ncol(x)),
            fn = neglogL_k, gr = neglogL_k_gradient, method = "BFGS",
            x = x, time = time, status = status, k = 1)


foo2 <- cmprsk::crr(time,status,x)

foo2$coef
foo$par
#neglogL_k_gradient(x, c(-0.1, 0.1, 0.1), time, status, 1)
# TODO: slow, it's faster to compute G for all times and then just select from that set iteratively

