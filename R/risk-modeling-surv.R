

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
  basesurv <- function(time.eval){
    
    stopifnot(length(time.eval) == 1)
    hdnom::glmnet_basesurv(time = time, event = status, lp = lp, 
                           times.eval = time.eval, centered = center)
    
  } # FUN
  
  # returns survival probability
  surv <- function(time.eval){
    
    basesurv_point <- basesurv(time.eval = time.eval)
    exp( -exp(lp) * basesurv_point$cumulative_base_hazard)
    
  } # FUN
  
  
  return(list(basesurv = basesurv, surv = surv))
  
} # FUN

# TODO: In baseline.risk of risk.modeling, there is error in coef: no var.kept there when naming the variables! The @x is wrong, as only nonzero gets printed!
# TODO: email Trevor about cryptic error in predict and coef.obj@i+1
# TODO: center argument
# TODO: return basesurv everywhere, put Breslow in computation of non-penalized survival!
# TODO: X needs column names (also in non-cox functions)
# TODO: the fact that response return is a function will cause issues in get.benefit()