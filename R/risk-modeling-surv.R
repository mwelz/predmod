
#' baseline survival w/o competing risks, with penalty
#' 
#' @param X covariate matrix
#' @param status binary status. 1 if failure, 0 if survived
#' @param time right-censored time at risk
#' @param alpha Alpha in glmnet
#' @param center Shall baseline survival be centered? Default is \code{FALSE}.
#' @param ... Additional parameters to be passed to \code{\link[glmnet]{cv.glmnet}}
#' 
#' @noRd
baseline.risk.surv.NoCmprsk.penalized <- function(X, status, time, 
                                                  alpha = 1, center = FALSE,
                                                  ...){
  
  # prepare dependent variable
  y <- survival::Surv(time = time, event = status)
  
  # fit models
  model.obj <- glmnet::cv.glmnet(x = X, y = y,
                                 family = "cox",
                                 type.measure = "deviance",
                                 alpha = alpha, ...) 
  
  # get linear predictor at best lambda: x'b
  #lp <- as.numeric(glmnet::predict.glmnet(
  #  model.obj,
  #  newx = data.matrix(X),
  #  s = "lambda.min", type = "link"
  #)) # throws cryptic error, so workaround required
  
  coef.obj <- glmnet::coef.glmnet(model.obj, s = "lambda.min")
  lp <- as.numeric(X[,coef.obj@i+1, drop = FALSE] %*% coef.obj@x) 
  
  # estimate survival 
  surv <- survival(time = time, status = status, lp = lp, center = center)
  
  # get coefficients
  coefs.obj        <- glmnet::coef.glmnet(model.obj, s = "lambda.min")
  kept.vars        <- coefs.obj@i
  coefs            <- as.matrix(coefs.obj)
  colnames(coefs)  <- "Estimate"
  
  return(list(
    model.obj = model.obj,
    lambda.min = model.obj$lambda,
    linear.predictor = lp,
    response = surv$surv,
    coefficients = coefs,
    retained.variables = colnames(X)[kept.vars],
    basesurv = surv$basesurv
  ))
  
} # FUN



#' baseline survival w/o competing risks, no penalty
#' 
#' @param X covariate matrix
#' @param status binary status. 1 if failure, 0 if survived
#' @param time right-censored time at risk
#' 
#' @noRd
baseline.risk.surv.NoCmprsk.ordinary <- function(X, status, time){
  
  # fit Cox PH model
  model.obj <- survival::coxph(survival::Surv(time = time, event = status)~.,
                               data = data.frame(time, status, X))
  
  # get linear predictor
  lp <- as.numeric(X %*% model.obj$coefficients)
  
  # survival probabilities is response. A function here
  f <- function(time){
    stopifnot(length(time) == 1) 
    as.numeric(summary(survival::survfit(model.obj, newdata = data.frame(X)), time = time)$surv)
  } # FUN
  
  # coefficients
  coefs            <- matrix(model.obj$coefficients)
  rownames(coefs)  <- colnames(X)
  colnames(coefs)  <- "Estimate"
  
  # return
  return(list(
    model.obj = model.obj,
    lambda.min = NULL,
    linear.predictor = lp,
    response = f,
    coefficients = coefs,
    retained.variables = colnames(X)
  ))
  
} # FUN


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