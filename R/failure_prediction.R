#' Function to predict the exact failure times
#' @param x A predmod object. TODO
#' @param time A vector of times along which to compute the survival curves
#'
#' @export
failure <- function(x, time = NULL)
{
  stopifnot(inherits(x, what = "risk_model_surv"))
  
  if(is.null(time)){
    time. <- x$inputs$time
  } else{
    time. <- time
  } # IF
  
  # extract estimated survival functions
  S_reg <- x$funs$regular$surv
  S_rev <- x$funs$counterfactual$surv
  
  # get sorted unique failure times to obtain survival curves
  time_unique    <- sort(unique(time.), decreasing = FALSE)
  surv_curve_reg <- S_reg(time_unique)
  surv_curve_rev <- S_rev(time_unique)
  
  # use survival curves to estimate potential failure times
  fail_reg <- expected_survival(S.hat = surv_curve_reg, Y.grid = time_unique)
  fail_rev <- expected_survival(S.hat = surv_curve_rev, Y.grid = time_unique)
  
  # calculate predicted benefits
  benefits <- get_predicted_benefits(risk_reg = fail_reg, 
                                     risk_rev = fail_rev,
                                     w = x$inputs$w)
  
  return(list(
    regular = fail_reg,
    counterfactual = fail_rev,
    benefits = benefits
  ))
  
} # FUN
