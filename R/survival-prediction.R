

#' estimates survival by using Breslow estimator
#' 
#' @param time A vector of times that was used to fit a Cox PH model
#' @param status A vector of mortality status that was used to fit a Cox PH model
#' @param lp A linear predictor, obtained from a Cox PH model
#' @param center Shall baseline survival be centered? Default is \code{FALSE}.
#' 
#' @noRd
survival <- function(time, status, lp, center = FALSE){
  
  # Breslow baseline survival function 
  basesurv <- function(time_eval){
    
    obj <- hdnom::glmnet_basesurv(
      time = time, event = status, lp = lp, 
      times.eval = time_eval, centered = center)
    obj$base_surv
    
  } # FUN
  
  # returns survival probability
  surv <- function(time_eval){
    
    t(outer(basesurv(time_eval = time_eval), exp(lp), FUN = "^"))
    
  } # FUN
  
  
  return(list(basesurv = basesurv, surv = surv))
  
} # FUN


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


#' creates helper object for predictions in Fine-Gray models
#' 
#' @param time time at risk
#' @param status Status. Zero if survived
#' @param k Status of interest. Default is 1.
#' @return A list
#' 
#' @export
prep_predict <- function(time, status, k = 1){
  
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
      
      # return
      list(Rk = Rk, zw = zw)
      
    }) # LAPPLY
  
  # naming
  names(Rk_list) <- paste0("idx", Ik)
  
  # return
  return(list(Ik = Ik, Rk_list = Rk_list, G = G, KM = km))
  
} # FUN


# estimates Fine-Gray Breslow baseline hazard function
basehaz_cmprsk <- function(time, status, lp, time_eval, prep_predict_object = NULL, failcode = 1){
  
  # prepare prediction
  if(is.null(prep_predict_object)){
    prep <- prep_predict(time = time, status = status, k = failcode)
  } else{
    prep <- prep_predict_object
  }
  
  # get Rk_list and Ik
  Ik <- prep$Ik
  Rk_list <- prep$Rk_list
  
  # get w_i(T_i)*N_{i,k}(t) for all i
  wn_mat <- 1L * outer(time, time_eval, FUN = "<=")
  
  # get exp(lp) on all samples
  explp <- exp(as.numeric(lp))
  
  # predict baseline survival
  mat <- matrix(sapply(Ik, function(i){
    
    zw  <- Rk_list[[paste0("idx",i)]]$zw
    Rk  <- Rk_list[[paste0("idx",i)]]$Rk
    wn_mat[i,,drop = FALSE] / sum(zw * explp[Rk])
    
  }), nrow = length(time_eval)) # SAPPLY
  
  rowSums(mat)
  
} # FUN

#' estimates survival in subdistributional hazard models
#' 
#' @param time A vector of times that was used to fit a cmprsk model
#' @param status A vector of mortality status that was used to fit a cmprsk model
#' @param lp A linear predictor, obtained from a cmprsk model
#' @param failcode The failure typr of interest that was used to fit a cmprsk model
#' @param prep_predict_object Output of \code{\link{prep_predict}}.
#' 
#' @noRd
survival_cmprsk <- function(time, status, lp, prep_predict_object = NULL, failcode = 1){
  
  # Fine-Gray Breslow baseline survival function
  basesurv <- function(time_eval){
    
    H0k <- basehaz_cmprsk(time = time, status = status, lp = lp, 
                          time_eval = time_eval, 
                          prep_predict_object = prep_predict_object, failcode = failcode)
    exp(-H0k) 
    
  } # FUN
  
  # function that returns survival curve of cause k
  surv <- function(time_eval){
    
    t(outer(basesurv(time_eval = time_eval), exp(lp), FUN = "^"))

  } # FUN
  
  # return
  return(list(basesurv = basesurv, surv = surv))
  
} # FUN


# TODO: In baseline.risk of risk.modeling, there is error in coef: no var.kept there when naming the variables! The @x is wrong, as only nonzero gets printed!
# TODO: email Trevor about cryptic error in predict and coef.obj@i+1
# TODO: center argument
# TODO: return basesurv everywhere, put Breslow in computation of non-penalized survival!
# TODO: X needs column names (also in non-cox functions)
# TODO: the fact that response return is a function will cause issues in get.benefit()