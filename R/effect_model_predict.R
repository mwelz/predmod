#' Predict method for a \code{effect_model_crss} object
#' 
#' @param object A \code{effect_model_crss} object.
#' @param newX A numeric matrix at which predictions should be performed
#' @param neww Treatment assignment
#' @param ... Additional parameters to be passed down
#' 
#' @return A matrix of risk predictions
#' 
#' @export
predict.effect_model <- function(object,
                                 neww, 
                                 newX,
                                 ...)
{
  
  ################################### all of the below is equal to predict.risk_model(),
  #consider merging!
  
  ## input checks
  if(!inherits(x = object, what = "effect_model_crss", which = FALSE))
  {
    stop("object must be an instance of effect_model_crss()")
  }
  
  ## check correctness of inpurs
  InputChecks_W(neww)
  InputChecks_newX(newX)
  InputChecks_equal.length2(neww, newX)
  newX <- check_and_adjust_newX(newX = newX, 
                                covariates = colnames(object$inputs$X))
   
  ## predict
  predict_effect_model_NoChecks(object = object,
                                neww = neww, 
                                newX = newX,
                                 ... = ...)
  
} # FUN



predict_effect_model_NoChecks <- function(object,
                                          neww, 
                                          newX,
                                          ...)
{
  
  
  ## get interacted variables
  interacted <- object$inputs$interacted
  
  if(is.null(interacted)){
    
    ## no variables are interacted
    
    # get full matrix
    X_full <- cbind(newX, w = neww)
    
    # prepare regressor matrix with reversed w
    X_full_rev <- cbind(newX, w = ifelse(neww == 1, 0, 1))
    
  } else{
    
    ## at least one variable is interacted
    
    # get full matrix
    X_full <- interacted_matrix(X = newX, w = neww, interacted = interacted)
    
    # prepare regressor matrix with reversed w
    X_full_rev <- interacted_matrix(X = newX, w = ifelse(neww == 1, 0, 1), interacted = interacted)
    
  } # IF
  
  
  ## get the retained variables
  retained <- rownames(object$coefficients$reduced)
  
  if("(Intercept)" %in% retained)
  {
    retained0 <- retained[-1L]
    intercept <- object$coefficients$reduced["(Intercept)", "Estimate"]
  } else{
    intercept <- 0.0
  }
  
  if(identical(length(retained0), 0L))
  {
    # case 1: no variables were retained
    risk_reg <- risk_rev <- rep(intercept, nrow(newX))
  } else
  {
    # case 2: at least one variable was retained
    X_full0 <- X_full[,retained0, drop = FALSE]
    X_full_rev0 <- X_full_rev[,retained0,drop = FALSE]
    
    colnames(X_full0) <- colnames(X_full_rev0) <- retained0
    
    ## predict risk with regular w... 
    risk_reg <- as.numeric(
      stats::plogis(cbind(1.0, X_full0) %*% object$models$reduced$coefficients)
    )
    
    ## ... and with flipped w
    risk_rev <- as.numeric(
      stats::plogis(cbind(1.0, X_full_rev0) %*% object$models$reduced$coefficients)
    )
    
  } # IF
  
  # calculate predicted benefits
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