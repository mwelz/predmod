#' Predict method for a \code{predmod_crss} object
#' 
#' @param object A \code{predmod_crss} object.
#' @param neww Treatment assignment
#' @param newz TODO
#' @param ... Additional parameters to be passed down
#' 
#' @return A matrix of risk predictions
#' 
#' @export
predict.risk_model <- function(object, neww, newz, ...)
{
  ## input checks
  if(!inherits(x = object, what = "risk_model_crss", which = FALSE))
  {
    stop("object must be an instance of risk_model_crss()")
  }
  
  ## check correctness of neww and newz
  InputChecks_W(neww)
  stopifnot(is.numeric(newz))
  InputChecks_equal.length2(neww, newz)
  
  ## predict
  predict_risk_model_NoChecks(object = object,
                              neww = neww, 
                              newz = newz, ... = ...)
  
} # FUN


predict_risk_model_NoChecks <- function(object, neww, newz, ...)
{
  ## get accepted model object from 2nd stage
  mod <- object$models$stage2$model$accepted
  
  ## flip w
  w_flipped <- ifelse(neww == 1, 0, 1)
  
  ## prepare X for second stage
  if(object$models$stage2$decision == "reduced"){
    
    # prepare X for second stage...
    X_stage2 <- cbind(w = neww, z = as.numeric(newz))
    
    # ... and its counterpart with flipped w
    X_stage2_rev <- cbind(w = w_flipped, z = as.numeric(newz))
    
  } else{
    
    # prepare X for second stage...
    X_stage2 <- cbind(w = neww,
                      z = as.numeric(newz), 
                      w.z = as.numeric(neww * newz))
    
    # ... and its counterpart with flipped w
    X_stage2_rev <- cbind(w = w_flipped,
                          z = as.numeric(newz),
                          w.z = as.numeric(w_flipped * newz))
    
  } # IF
  
  
  ## predict risk with regular w... 
  risk_reg <- as.numeric(
    stats::plogis(cbind(1.0, X_stage2) %*% mod$coefficients)
  )
  
  ## ... and with flipped w
  risk_rev <- as.numeric(
    stats::plogis(cbind(1.0, X_stage2_rev) %*% mod$coefficients)
  )
  
  
  ## calculate predicted benefits
  benefits <- get_predicted_benefits_NoChecks(
    risk_reg = risk_reg, 
    risk_rev = risk_rev,
    w = neww)
  
  return(cbind(
    benefit_absolute = benefits$absolute,
    benefit_relative = benefits$relative,
    risk_regular = risk_reg,
    risk_counterfactual = risk_rev
  ))
} # FUN
