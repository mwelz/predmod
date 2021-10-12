rm(list = ls()) ; cat("\014")

km_fit <- function(time, status){
  
  # status should be binary
  survival::survfit(survival::Surv(time, status) ~ 1, 
                    data = data.frame(time, status),
                    type = "kaplan-meier")
} # FUN


km_predict_survival <- function(x, time){
  
  # x is returned by km_fit()
  summary(x, time = time)$surv
  
} # FUN


optim_prep_pshm <- function(x, time, status, k){
  
  # get Kaplan-Meier (KM) estimator of censoring time by using {time, 1 - delta}
  km <- km_fit(time = time, status = 1*(status == 0))
  
  # evaluate KM estimator at all times
  G <- sapply(1:length(time), function(i) {
    km_predict_survival(x = km, time = time[i])
    })
  
  # get index set Ik
  Ik <- which(status == k)
  
  # get Rk_list
  Rk_list <- 
    lapply(Ik, function(i){
      
      # get set R_{i,k}
      Rk <- which(time >= time[i] | (0 < status & status != k))
      
      # get Z_j(Y_i) * w_j(Y_i) for all j \in R_{i,k}
      zw <- sapply(Rk, function(j){
        ifelse(test = (time[j] < time[i]) & (0 < status[j]) & (status[j] != k), 
               yes = G[i] / G[j], no = 1.0)
      })
      
      # get Z_j(Y_i) * w_j(Y_i) * X_j for all j \in R_{i,k}
      zwx <- x[Rk,] * zw
      
      # return
      list(Rk = Rk, zw = zw, zwx = zwx)
      
    }) # LAPPLY
  
  # naming
  names(Rk_list) <- paste0("idx", Ik)
  
  # return
  return(list(Ik = Ik, Rk_list = Rk_list, G = G, KM = km))
  
} # FUN


neglogL <- function(x, beta, optim_prep_pshm_object){
  
  # get X*beta and exp() thereof
  xb  <- as.numeric(x %*% beta)
  exb <- exp(xb)
  
  # get Rk_list and Ik
  Ik <- optim_prep_pshm_object$Ik
  Rk_list <- optim_prep_pshm_object$Rk_list
  
  # get the log likelihood
  -sum(sapply(Ik, function(i){
    
    Rk  <- Rk_list[[paste0("idx",i)]]$Rk
    zw  <- Rk_list[[paste0("idx",i)]]$zw

    xb[i] - log(sum(zw * exb[Rk]))
    
  })) # SAPPLY
  
} # FUN


neglogL_gradient <- function(x, beta, optim_prep_pshm_object){
  
  # get X*beta and exp() thereof
  xb  <- as.numeric(x %*% beta)
  exb <- exp(xb)
  p   <- ncol(x)
  
  # get Rk_list and Ik
  Ik <- optim_prep_pshm_object$Ik
  Rk_list <- optim_prep_pshm_object$Rk_list
  
  # evaluate the gradient
  -rowSums(
    
    sapply(Ik, function(i){
    
    Rk  <- Rk_list[[paste0("idx",i)]]$Rk
    zw  <- Rk_list[[paste0("idx",i)]]$zw
    zwx <- Rk_list[[paste0("idx",i)]]$zwx
    
    x[i,] - colSums(matrix(zwx * exb[Rk], ncol = p)) / (sum(zw * exb[Rk]))
    
  })) # SAPPLY
  
} # FUN


#' fitting a proportional subdistribution hazard model as in Fine and Gray (1999, JASA)
#' 
#' @param x Matrix of covariates.
#' @param time Vector of failure/censoring times
#' @param status Vector with a unique code (nonzero) for each failure type and a separate code (equality to zero) for censored observations
#' @param k fit for k-th failure type.
#' 
#' TODO: covariance estimate, dealing with duplicate failure times, Breslow-type estimator for prediction
#' 
#' @export
pshm <- function(x, time, status, k = 1){
  
  x <- as.matrix(x)
  if(is.null(colnames(x))) colnames(x) <- paste0("V", 1:ncol(x))
  
  # prepare and perform optimization routine
  optim_prep_pshm_object <- optim_prep_pshm(x = x, time = time, status = status, k = k)
  opt <- optim(par = rep(0.0, ncol(x)),
               fn = neglogL,
               gr = neglogL_gradient, 
               method = "BFGS",
               x = x, 
               optim_prep_pshm_object = optim_prep_pshm_object)
  
  return(structure(list(
    coef = structure(opt$par, names = colnames(x)),
    loglik = -opt$value,
    km = optim_prep_pshm_object$KM
  ), class = "pshm"))
} # FUN


########################################

# simulated data to test
set.seed(10)
time <- rexp(200)
status <- sample(0:2,200,replace=TRUE)
x <- matrix(runif(600),nrow=200)
dimnames(x)[[2]] <- c('x1','x2','x3')

pshm.obj <- pshm(x = x, time = time, status = status, k = 1)
pshm.obj$coef

# for comparison:
cmprsk::crr(ftime = time, fstatus = status, cov1 = x)$coef # accurate
