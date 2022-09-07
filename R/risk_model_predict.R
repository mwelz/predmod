#' Predict method for a \code{predmod_crss} object
#' 
#' @param object A \code{predmod_crss} object.
#' @param newX A numeric matrix at which predictions should be performed
#' @param neww Treatment assignment
#' @param newz TODO
#' @param ... Additional parameters to be passed down
#' 
#' @return A matrix of risk predictions
#' 
#' @export
predict.risk_model <- function(object, neww, newz = NULL, newX = NULL, ...)
{
  ## input checks
  if(!inherits(x = object, what = "risk_model_crss", which = FALSE))
  {
    stop("object must be an instance of risk_model_crss()")
  }
  
  ## check correctness of neww
  InputChecks_W(neww)
  
  ## if no newz is passed, we need both newX and baseline model
  if(is.null(newz))
  {
    if(is.null(object$models$baseline))
    {
      ## case 1: no baseline model => stop
      stop(paste0("No baseline risk model was fitted in 'object'. ",
                  " Therefore, please provide a vector 'newz'"), 
           call. = FALSE)
    } else if(is.null(newX)){
      
      ## case 2: no newX => stop
      stop("Please pass a matrix 'newX' to generate baseline risk predictions", 
           call. = FALSE) 
    } else{
      
      ## case 3: both baseline model and newX passed
      # check if newX is correctly specified
      InputChecks_newX(newX)
      InputChecks_equal.length2(neww, newX)
      newX <- check_and_adjust_newX(newX = newX, 
                                    covariates = object$inputs$covariates)
      #InputChecks_newX_X(newX = newX,
      #                   object = object$models$baseline$, 
      #                   survival = FALSE)
    } # IF
  } # IF
  
  ## predict
  predict_risk_model_NoChecks(object = object,
                              neww = neww, 
                              newX = newX,
                              newz = newz, ... = ...)
  
} # FUN


predict_risk_model_NoChecks <- function(object, neww, newX, newz, ...)
{
  ## get accepted model object from 2nd stage
  mod <- object$models$stage2$model$accepted
  
  ## if no newz provided, get it via the baseline risk model
  if(is.null(newz))
  {
    
    ## get baseline risk
    risk0 <- predict_baseline_crss_NoChecks(
      object = object$models$baseline, 
      newX = newX, ... = ...)
    
    ## take as z the linear predictor of risk
    z <- stats::qlogis(p = risk0)
    
  } else{
    z <- newz
  } # IF
  
  ## flip w
  w_flipped <- ifelse(neww == 1, 0, 1)
  
  ## prepare X for second stage
  if(object$models$stage2$decision == "reduced"){
    
    # prepare X for second stage...
    X_stage2 <- cbind(w = neww, z = as.numeric(z))
    
    # ... and its counterpart with flipped w
    X_stage2_rev <- cbind(w = w_flipped, z = as.numeric(z))
    
  } else{
    
    # prepare X for second stage...
    X_stage2 <- cbind(w = neww,
                      z = as.numeric(z), 
                      w.z = as.numeric(neww * z))
    
    # ... and its counterpart with flipped w
    X_stage2_rev <- cbind(w = w_flipped,
                          z = as.numeric(z),
                          w.z = as.numeric(w_flipped * z))
    
  } # IF
  
  
  ## predict risk with regular w... 
  risk_reg <- as.numeric(
    stats::predict.glm(object = mod, 
                       newdata = as.data.frame(X_stage2),
                       type = "response"))
  
  ## ... and with flipped w
  risk_rev <- as.numeric(
    stats::predict.glm(object = mod, 
                       newdata = as.data.frame(X_stage2_rev),
                       type = "response"))
  
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



#' Predict method for a \code{predmod_surv} object
#' 
#' @param object A \code{predmod_surv} object.
#' @param newX A numeric matrix at which predictions should be performed
#' @param neww Treatment assignment
#' @param newz TODO
#' @param time_eval The time at which baseline risk shall be predicted. Must be nonnegative numeric vector of length one.
#' @param ... Additional parameters to be passed down
#' 
#' @return A matrix of risk predictions
#' 
#' @export
predict.risk_model_survival <- function(object, neww, time_eval, newz = NULL, newX = NULL, ...)
{
  ## input checks
  if(!inherits(x = object, what = "predmod_surv", which = FALSE))
  {
    stop("object must be an instance of predmod_surv()")
  }
  
  ## check correctness of neww and time_eval
  InputChecks_W(neww)
  stopifnot(length(time_eval) == 1L & is.numeric(time_eval))
  
  ## if no newz is passed, we need both newX and baseline model
  if(is.null(newz))
  {
    if(is.null(object$models$baseline))
    {
      ## case 1: no baseline model => stop
      stop(paste0("No baseline risk model was fitted in 'object'. ",
                  " Therefore, please provide a vector 'newz'"), 
           call. = FALSE)
    } else if(is.null(newX)){
      
      ## case 2: no newX => stop
      stop("Please pass a matrix 'newX' to generate baseline risk predictions", 
           call. = FALSE) 
    } else{
      
      ## case 3: both baseline model and passed
      # check if newX is correctly specified
      InputChecks_newX(newX)
      InputChecks_equal.length2(neww, newX)
      InputChecks_newX_X(newX = newX,
                         object = object$models$baseline, 
                         survival = TRUE)
    } # IF
  } # IF
  
  ## predict
  predict_risk_survival_NoChecks(object    = object,
                                 neww      = neww, 
                                 time_eval = time_eval,
                                 newX      = newX,
                                 newz      = newz, ... = ...)
  
} # FUN


predict_risk_survival_NoChecks <- function(object, neww, time_eval, newX, newz, ...)
{
  ## get model object from 2nd stage
  mod <- object$models$stage2$model
  
  ## if no newz provided, get it via the baseline risk model
  if(is.null(newz))
  {
    
    ## get baseline risk
    risk0 <- predict_baseline_survival_NoChecks(
      object = object$models$baseline,
      time_eval = time_eval,
      newX = newX, ... = ...)
    
    ## take as z the linear predictor of risk
    z <- stats::qlogis(p = risk0)
    
  } else{
    z <- newz
  } # IF
  
  ## flip w
  w_flipped <- ifelse(neww == 1, 0, 1)
  
  ## prepare X for second stage
  if(object$inputs$constant){
    
    # prepare X for second stage...
    X_stage2 <- cbind(w = neww, z = as.numeric(z))
    
    # ... and its counterpart with flipped w
    X_stage2_rev <- cbind(w = w_flipped, z = as.numeric(z))
    
  } else{
    
    # prepare X for second stage...
    X_stage2 <- cbind(w = neww,
                      z = as.numeric(z), 
                      w.z = as.numeric(neww * z))
    
    # ... and its counterpart with flipped w
    X_stage2_rev <- cbind(w = w_flipped,
                          z = as.numeric(z),
                          w.z = as.numeric(w_flipped * z))
    
  } # IF
  
  stop("continue with the prediction method of survival risk model")
  
} # FUN
