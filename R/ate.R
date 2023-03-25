#' computes the average treatment effect (ATE) of a predmod object
#' 
#' @param x A predmod object
#' @param subset The indices of the subgroup of interest
#' @param relative Shall relative ATE be calculated?
#' 
#' @export
average_treatment_effect <- function(x, 
                                     subset = NULL, 
                                     relative = FALSE,
                                     neww = NULL, 
                                     newX = NULL,
                                     newz = NULL,
                                     shrunk = FALSE){
  
  stopifnot(inherits(x, what = c("risk_model_crss",
                                 "effect_model_crss",
                                 "causal_forest")))
  clss <- class(x)
  
  # prepare subset object
  if(!is.null(subset)){
    stopifnot(is.numeric(subset))
  } else {
    subset <- seq_along(x$inputs$w)
  } # IF
  
  nullX <- is.null(newX)
  nullw <- is.null(neww)
  nullz <- is.null(newz)
  
 if(clss == "effect_model_crss")
 {
   if(nullX && nullw)
   {
     ## no new data passed, so use data that model was fitted on
     X0 <- x$inputs$X
     w0 <- x$inputs$w
     
   } else if(!nullX && !nullw)
   {
     ## new data passed, so use the new data
     # checks for the correct format of the new data will be done in
     # predict.effect_model(), which will be called later
     w0 <- neww
     X0 <- newX
   } else{
     stop("newX and neww must either be both NULL or non-NULL for effect models")
   }
   
   # effect models don't have a 'z' component, so set it to NULL
   z0 <- NULL
   
 } else #if(clss == "risk_model_crss")
 {
   
   if(nullw && nullz)
   {
     ## no new data passed, so use data that model was fitted on
     # that is, recover z and w: these are always present in 2nd stage of risk model
     z0 <- x$inputs$z
     w0 <- x$inputs$w
   } else if(!nullw && !nullz)
   {
     ## new data passed, so use the new data
     # checks for the correct format of the new data will be done in
     # predict.risk_model(), which will be called later
     z0 <- newz
     w0 <- neww
   } else{
     stop("newz and neww must either be both NULL or non-NULL for risk models")
   }
   
   # risk models don't have a 'X0' component, so set it to NULL
   X0 <- NULL
    
 } # IF clss
  
  ## call the main function
  average_treatment_effect_crss_NoChecks(x = x, 
                                         subset = subset, 
                                         relative = relative,
                                         neww = w0, 
                                         newX = X0,
                                         newz = z0,
                                         shrunk = shrunk)
  
} # FUN

# all arguments are non-null, come from average_treatment_effect()
average_treatment_effect_crss_NoChecks <- function(x, 
                                                   subset , 
                                                   relative,
                                                   neww, 
                                                   newX,
                                                   newz,
                                                   shrunk)
{
  if(relative)
  {
    ATE_relative_crss(x = x, 
                      subset = subset, 
                      neww = neww, 
                      newX = newX,
                      newz = newz,
                      shrunk = shrunk)
  } else{
    ATE_absolute_crss(x = x, 
                      subset = subset, 
                      neww = neww, 
                      newX = newX,
                      newz = newz,
                      shrunk = shrunk)
  } # IF
} # FUN


# all arguments are non-null, come from average_treatment_effect()
ATE_absolute_crss <- function(x, 
                              subset, 
                              neww, 
                              newX,
                              newz,
                              shrunk)
{
  # note: we don't check newz or newX to be NULL because this would cause
  # incompatibilities in effect models (which don't have z) or risk
  # models (which don't have X)
  # The corresponding checks and error prompts are in  average_treatment_effect()
  if(inherits(x = x, what = "risk_model_crss"))
  {
    pred <- predict.risk_model(object = x, 
                               neww = neww, 
                               newz = newz)
  } else if(inherits(x = x, what = "effect_model_crss"))
  {
    pred <- predict.effect_model(object = x, 
                                 neww = neww, 
                                 newX = newX, 
                                 shrunk = shrunk)
  } else{
    stop("predict methods for GRF aren't yet implemented")
    pred <- predict.grf_model(object = x, 
                              newX = newX)
  } # IF inherits
  
  # get risks
  x_reg <- pred[,"risk_regular"]
  x_rev <- pred[,"risk_counterfactual"]
  w     <- neww
  
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



# all arguments are non-null, come from average_treatment_effect()
ATE_relative_crss <- function(x, 
                              subset, 
                              neww, 
                              newX,
                              newz,
                              shrunk)
{
  # note: we don't check newz or newX to be NULL because this would cause
  # incompatibilities in effect models (which don't have z) or risk
  # models (which don't have X)
  # The corresponding checks and error prompts are in  average_treatment_effect()
  if(inherits(x = x, what = "risk_model_crss"))
  {
    pred <- predict.risk_model(object = x, 
                               neww = neww, 
                               newz = newz)
  } else if(inherits(x = x, what = "effect_model_crss"))
  {
    pred <- predict.effect_model(object = x, 
                                 neww = neww, 
                                 newX = newX,
                                 shrunk = shrunk)
  } else{
    stop("predict methods for GRF aren't yet implemented")
    pred <- predict.grf_model(object = x, 
                              newX = newX)
  } # IF inherits
  
  # get ATE
  ate <- mean(pred[subset,"benefit_relative"])
  
  # no SE can be computed for relative risk TODO: maybe via bootstrap
  sderr <- NA_real_
  
  # return
  return(structure(c(ate, sderr), names = c("ATE", "Std. Error")))
  
} # FUN
