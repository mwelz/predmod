rm(list = ls()) ; cat("\014")

optim_prep_phm <- function(x, time, status){
  
  # get index set I
  I <- which(status == 1)
  
  # get R_list
  R_list <- 
    lapply(I, function(i){
      
      # get set R_{i,k}
      which(time >= time[i])
      
    }) # LAPPLY
  
  # naming
  names(R_list) <- paste0("idx", I)
  
  # return
  return(list(I = I, R_list = R_list))
  
} # FUN


neglogL <- function(x, beta, optim_prep_phm_object){
  
  # get X*beta and exp() thereof
  xb  <- as.numeric(x %*% beta)
  exb <- exp(xb)
  
  # get Rk_list and Ik
  I <- optim_prep_phm_object$I
  R_list <- optim_prep_phm_object$R_list
  
  # get the log likelihood
  -sum(sapply(I, function(i){
    
    R  <- R_list[[paste0("idx",i)]]
    
    xb[i] - log(sum(exb[R]))
    
  })) # SAPPLY
  
} # FUN




neglogL_gradient <- function(x, beta, optim_prep_phm_object){
  
  # get X*beta and exp() thereof
  xb  <- as.numeric(x %*% beta)
  exb <- exp(xb)
  p   <- ncol(x)
  
  # get Rk_list and Ik
  I <- optim_prep_phm_object$I
  R_list <- optim_prep_phm_object$R_list
  
  # get the gradient, evaluated at beta
  -rowSums(sapply(I, function(i){
    
    R  <- R_list[[paste0("idx",i)]]
    
    x[i,] - colSums(matrix(x[R,] * exb[R], ncol = p)) / sum(exb[R])
    
  })) # SAPPLY
  
} # FUN



#' fitting a proportional hazard model Cox (1972, JRSSB)
#' 
#' @param x Matrix of covariates.
#' @param time Vector of failure/censoring times
#' @param status Binary vector. Equal to to zero for censoring times.
#' 
#' TODO: covariance estimate, dealing with duplicate failure times, prediction
#' 
#' @export
phm <- function(x, time, status){
  
  x <- as.matrix(x)
  if(is.null(colnames(x))) colnames(x) <- paste0("V", 1:ncol(x))
  
  # prepare and perform optimization routine
  optim_prep_phm_object <- optim_prep_phm(x = x, time = time, status = status)
  opt <- optim(par = rep(0.0, ncol(x)),
               fn = neglogL,
               gr = neglogL_gradient, 
               method = "BFGS",
               x = x, 
               optim_prep_phm_object = optim_prep_phm_object)
  
  return(structure(list(
    coef = structure(opt$par, names = colnames(x)),
    loglik = -opt$value
  ), class = "phm"))
} # FUN


km_fit <- function(time, status, type = "kaplan-meier"){
  
  # status should be binary
  survival::survfit(survival::Surv(time, status) ~ 1, 
                    data = data.frame(time, status),
                    type = type)
} # FUN


km_predict_survival <- function(x, time){
  
  # x is returned by km_fit()
  summary(x, time = time)$surv
  
} # FUN


############################

lung <- survival::lung

# 1 means having the event
lung$status <- ifelse(lung$status == 1, 1, 0)

# for now, kick out duplicate times
lung <- lung[!duplicated(lung$time),]

x <- cbind(sex = lung$sex, age = lung$age)
time <- lung$time
status <- lung$status

phm(x = x, time = time, status = status)$coef

# survival predictions at time[1]
phm_fit <- phm(x = x, time = time, status = status)
km <- km_fit(time, status, type = "kaplan-meier")
baseline_survival <- km_predict_survival(km, time[1])
betahat <- phm_fit$coef

pred <- baseline_survival^exp(as.numeric(x %*% betahat))



# comparison
library(survival)
check <- coxph(Surv(time, status) ~ sex + age, data = lung )
check_pred <- summary(survfit(check, newdata = data.frame(x)), time=time[1])
as.numeric(check_pred$surv)

bhaz=basehaz(check, centered = T)
exp(-bhaz[bhaz$time == time[1],1])

baseline_survival # we get close, but not exactly the same :(

