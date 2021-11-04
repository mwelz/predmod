# baseline risk prediction via penalized logistic regression (treatment assignment is not used)
# 
# @param X matrix of data frame of covariates
# @param y vector of binary responses
# @param alpha the alpha as in glmnet. Defauly is 1 (= Lasso). If NULL, no penalty is applied.
# 
# @export
# baseline.risk <- function(X, y, alpha = 1){
#   
#   if(is.null(alpha)){
#     
#     # fit the model
#     model.obj <- stats::glm(y~., 
#                             family = stats::binomial(link = "logit"), 
#                             data = data.frame(y, X))
#     
#     # get liner predictor etc.
#     lp        <- as.numeric(cbind(1,X) %*% model.obj$coefficients)
#     kept.vars <- 1:ncol(X)
#     
#     coefs            <- matrix(model.obj$coefficients)
#     rownames(coefs)  <- c("(Intercept)", colnames(X)[kept.vars])
#     colnames(coefs)  <- "Estimate"
#     lambda           <- NULL
#     
#   } else{
#     
#     X             <- as.matrix(X)
#     if(is.null(colnames(X))) colnames(X) <- paste0("X", 1:ncol(X))
#     model.obj    <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = alpha)
#     lambda        <- model.obj$lambda.min # minimizing lambda (needs to be fixed in second stage)
#     coefs.obj     <- glmnet::coef.glmnet(model.obj, s = "lambda.min")
#     kept.vars     <- coefs.obj@i
#     if(0 %in% kept.vars) kept.vars <- kept.vars[-which(kept.vars == 0)] # intercept is not a variable
#     coefs            <- as.matrix(coefs.obj@x)
#     rownames(coefs)  <- c("(Intercept)", colnames(X)[kept.vars])
#     colnames(coefs)  <- "Estimate"
#     X.retained       <- cbind(intercept = 1, X[,kept.vars,drop = FALSE])
#     lp               <- as.numeric(X.retained %*% coefs) # linear predictor
#     
#   } # IF
#   
#   
#   
#   return(list(
#     model.obj = model.obj,
#     lambda.min = lambda,
#     linear.predictor = lp,
#     response = stats::plogis(lp),
#     coefficients = coefs,
#     retained.variables = colnames(X)[kept.vars]
#   ))
# } # FUN


# baseline survival w/o competing risks, with penalty
# 
# @param X covariate matrix
# @param status binary status. 1 if failure, 0 if survived
# @param time right-censored time at risk
# @param alpha Alpha in glmnet
# @param center Shall baseline survival be centered? Default is \code{FALSE}.
# @param ... Additional parameters to be passed to \code{\link[glmnet]{cv.glmnet}}
# 
# @noRd
# baseline.risk.surv.NoCmprsk.penalized <- function(X, status, time, 
#                                                   alpha = 1, center = FALSE,
#                                                   ...){
#   
#   # prepare dependent variable
#   y <- survival::Surv(time = time, event = status)
#   
#   # fit models
#   model.obj <- glmnet::cv.glmnet(x = X, y = y,
#                                  family = "cox",
#                                  type.measure = "deviance",
#                                  alpha = alpha, ...) 
#   
#   # get linear predictor at best lambda: x'b
#   #lp <- as.numeric(glmnet::predict.glmnet(
#   #  model.obj,
#   #  newx = data.matrix(X),
#   #  s = "lambda.min", type = "link"
#   #)) # throws cryptic error, so workaround required
#   
#   coef.obj <- glmnet::coef.glmnet(model.obj, s = "lambda.min")
#   lp <- as.numeric(X[,coef.obj@i+1, drop = FALSE] %*% coef.obj@x) 
#   
#   # estimate survival 
#   surv <- survival(time = time, status = status, lp = lp, center = center)
#   
#   # get coefficients
#   coefs.obj        <- glmnet::coef.glmnet(model.obj, s = "lambda.min")
#   kept.vars        <- coefs.obj@i
#   coefs            <- as.matrix(coefs.obj)
#   colnames(coefs)  <- "Estimate"
#   
#   return(list(
#     model.obj = model.obj,
#     lambda.min = model.obj$lambda,
#     linear.predictor = lp,
#     response = surv$surv,
#     coefficients = coefs,
#     retained.variables = colnames(X)[kept.vars],
#     basesurv = surv$basesurv
#   ))
#   
# } # FUN



# baseline survival w/o competing risks, no penalty
# 
# @param X covariate matrix
# @param status binary status. 1 if failure, 0 if survived
# @param time right-censored time at risk
# @param center Shall baseline survival be centered? Default is \code{FALSE}.
# @param ... Additional parameters to be passed to \code{\link[glmnet]{cv.glmnet}}
# 
# @noRd
# baseline.risk.surv.NoCmprsk.ordinary <- function(X, status, time, center, ...){
#   
#   # fit Cox PH model
#   model.obj <- survival::coxph(survival::Surv(time = time, event = status)~.,
#                                data = data.frame(time, status, X), ...)
#   
#   # get linear predictor
#   lp <- as.numeric(X %*% model.obj$coefficients)
#   
#   # estimate survival with Breslow baseline. A function here 
#   surv <- survival(time = time, status = status, lp = lp, center = center)$surv
#   #as.numeric(summary(survival::survfit(model.obj, newdata = data.frame(X)), time = time)$surv)
#   
# 
#   # coefficients
#   coefs            <- matrix(model.obj$coefficients)
#   rownames(coefs)  <- colnames(X)
#   colnames(coefs)  <- "Estimate"
#   
#   # return
#   return(list(
#     model.obj = model.obj,
#     lambda.min = NULL,
#     linear.predictor = lp,
#     response = surv,
#     coefficients = coefs,
#     retained.variables = colnames(X)
#   ))
#   
# } # FUN



# # baseline survival w/ competing risks (Fine-Gray), no penalty
# # 
# # @param X covariate matrix
# # @param status binary status. 1 if failure, 0 if survived
# # @param time right-censored time at risk
# # @param failcode k
# # @param ... Additional parameters to be passed to \code{\link[glmnet]{cv.glmnet}}
# # 
# # @noRd
# baseline.risk.surv.cmprsk.ordinary <- function(X, status, time, failcode = 1, ...){
#   
#   # fit Cox PH model
#   model.obj <- cmprsk::crr(ftime = time, fstatus = status, 
#                            cov1 = X, failcode = failcode, cencode = 0)
#   
#   # get linear predictor
#   lp <- as.numeric(X %*% model.obj$coef)
#   
#   # estimate survival. returns a function
#   surv <- survival_cmprsk(time = time, status = status, lp = lp,
#                           prep_predict_object = NULL, k = failcode)$surv
#   
#   # coefficients
#   coefs            <- matrix(model.obj$coef)
#   rownames(coefs)  <- colnames(X)
#   colnames(coefs)  <- "Estimate"
#   
#   # return
#   return(list(
#     model.obj = model.obj,
#     lambda.min = NULL,
#     linear.predictor = lp,
#     response = surv,
#     coefficients = coefs,
#     retained.variables = colnames(X)
#   ))
#   
# } # FUN


# # baseline survival 
# # 
# # @param X covariate matrix
# # @param status binary status. 1 if failure, 0 if survived
# # @param time right-censored time at risk
# # @param alpha penalty
# # 
# # @export
# baseline.risk.survival <- function(X, status, time, alpha = NULL){
#   
#   if(length(unique(status)) > 2){
#     # competing risks
#     
#     if(is.null(alpha)){
#       
#       baseline.risk.surv.cmprsk.ordinary(X = X, status = status, time = time, failcode = 1)
#       
#     } else stop("Penalized cmprsk not yet implemented")
#     
#   } else{
#     # no competing risks
#     
#     if(is.null(alpha)){
#       
#       baseline.risk.surv.NoCmprsk.ordinary(X = X, status = status, time = time, center = FALSE)
#       
#     } else{
#       
#       baseline.risk.surv.NoCmprsk.penalized(X = X, status = status, 
#                                             time = time, alpha = alpha, center = FALSE)
#       
#     } # IF
#     
#   } # IF
#   
# } # FUN