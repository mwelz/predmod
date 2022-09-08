#' computes the average treatment effect (ATE) of a predmod object
#' 
#' @param x A predmod object
#' @param subset The indices of the subgroup of interest
#' @param relative Shall relative ATE be calculated?
#' @param time_eval Time at which to evaluate the failure risk predictions.
#' 
#' @export
average_treatment_effect <- function(x, 
                                     subset = NULL, 
                                     relative = FALSE,
                                     time_eval = NULL,
                                     neww = NULL, 
                                     newX = NULL,
                                     newz = NULL){
  
  stopifnot(inherits(x, what = c("risk_model_crss", "risk_model_surv",
                                 "effect_model_crss", "effect_model_surv",
                                 "grf_crss")))
  
  # prepare subset object
  if(!is.null(subset)){
    stopifnot(is.numeric(subset))
  } else {
    subset <- seq_along(x$inputs$w)
  } # IF
  
  
  ## another input check TODO: come up with better check
  #if(any(c(!is.null(neww), !is.null(newz), !is.null(newX), !is.null(subset))))
  #{
  #  stop("if one of neww, newX, newz, subset is not NULL, all must be non-null")
  #} # IF
  
  # call the correct function
  if(inherits(x, 
              what = c("risk_model_crss", "effect_model_crss", 
                       "grf_crss")))
  {
    average_treatment_effect_crss_NoChecks(x = x, 
                                           subset = subset, 
                                           relative = relative,
                                           neww = neww, 
                                           newX = newX,
                                           newz = newz)
  } else{
    stop("ATE for survival models not yet implemented")
    ## the below works, but is incomplete for external validation
    # average_treatment_effect_surv(x, 
    #                               subset = subset, 
    #                               relative = relative,
    #                               time_eval = time_eval)
  } # IF
} # FUN


average_treatment_effect_crss_NoChecks <- function(x, 
                                                   subset = NULL, 
                                                   relative = FALSE,
                                                   neww = NULL, 
                                                   newX = NULL,
                                                   newz = NULL)
{
  if(relative)
  {
    ATE_relative_crss(x = x, 
                      subset = subset, 
                      neww = neww, 
                      newX = newX,
                      newz = newz)
  } else{
    ATE_absolute_crss(x = x, 
                      subset = subset, 
                      neww = neww, 
                      newX = newX,
                      newz = newz)
  } # IF
} # FUN


ATE_absolute_crss <- function(x, 
                              subset = NULL, 
                              neww = NULL, 
                              newX = NULL,
                              newz = NULL)
{
  
  if(is.null(neww) && is.null(newX) && is.null(newz))
  {
    x_reg <- x$risk$regular
    x_rev <- x$risk$counterfactual
    w     <- x$inputs$w
  } else
  {
    
    if(inherits(x = x, what = "risk_model_crss"))
    {
      pred <- predict.risk_model(object = x, 
                                 neww = neww, 
                                 newz = newz, 
                                 newX = newX)
    } else if(inherits(x = x, what = "effect_model_crss"))
    {
      pred <- predict.effect_model(object = x, 
                                   neww = neww, 
                                   newX = newX)
    } else{
      pred <- predict.grf_model(object = x, 
                                newX = newX)
    } # IF inherits
    
    # get risks
    x_reg <- pred[,"risk_regular"]
    x_rev <- pred[,"risk_counterfactual"]
    w     <- neww
    
  } # IF is null
  
  # adjust signs to ensure that x_reg - x_rev = (predicted absolute benefit) 
  x_reg[w == 0] <- -x_reg[w == 0]
  x_rev[w == 0] <- -x_rev[w == 0]
  
  # perform t-test
  t <- stats::t.test(x = x_reg[subset], y = x_rev[subset], paired = FALSE, var.equal = FALSE)
  
  # get ATE
  ate <- unname(t$estimate[1] - t$estimate[2])
  sderr <- t$stderr
  
  # return
  return(structure(c(ate, sderr), names = c("ATE", "Std. Error")))
} # FUN




ATE_relative_crss <- function(x, 
                              subset = NULL, 
                              neww = NULL, 
                              newX = NULL,
                              newz = NULL)
{
  # relative effect: here we can simply use the direct estimated benefits
  if(is.null(neww) && is.null(newX) && is.null(newz))
  {
    ate <- mean(x$benefits$relative[subset])
  } else
  {
    
    if(inherits(x = x, what = "risk_model_crss"))
    {
      pred <- predict.risk_model(object = x, 
                                 neww = neww, 
                                 newz = newz, 
                                 newX = newX)
    } else if(inherits(x = x, what = "effect_model_crss"))
    {
      pred <- predict.effect_model(object = x, 
                                   neww = neww, 
                                   newX = newX)
    } else{
      pred <- predict.grf_model(object = x, 
                                newX = newX)
    } # IF inherits
    
    ate <- mean(pred[,"benefit_relative"])

  } # IF is null
  
  # no SE can be computed for relative risk TODO: maybe via bootstrap
  sderr <- NA_real_
  
  # return
  return(structure(c(ate, sderr), names = c("ATE", "Std. Error")))
  
} # FUN


## TODO: this is work in progress
average_treatment_effect_surv <- function(x, 
                                          subset = NULL, 
                                          relative = FALSE,
                                          time_eval = NULL)
{
  w <- x$inputs$w
  
  if(!relative)
  {
    # case 1: absolute effect: benefits concern risk at time of interest
    x_reg <- as.numeric(1.0 - x$funs$regular$surv(time_eval))
    x_rev <- as.numeric(1.0 - x$funs$counterfactual$surv(time_eval))
    
    # adjust signs to ensure that x_reg - x_rev = (predicted absolute benefit) 
    x_reg[w == 0] <- -x_reg[w == 0]
    x_rev[w == 0] <- -x_rev[w == 0]
    
    # perform t-test
    t <- stats::t.test(x = x_reg[subset], y = x_rev[subset], paired = FALSE, var.equal = FALSE)
    
    # get ATE
    ate <- unname(t$estimate[1] - t$estimate[2])
    sderr <- t$stderr
    
  } else{
    
    # case 2: relative effect: benefits concern risk at time of interest
    ate <- mean(x$benefits$relative[subset])
    
    # no SE can be computed for relative risk TODO: maybe via bootstrap
    sderr <- NA_real_
  } # IF
  
  # return
  return(structure(c(ate, sderr), names = c("ATE", "Std. Error")))
} # FUN
