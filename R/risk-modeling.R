#' baseline risk prediction via penalized logistic regression (treatment assignment is not used)
#' 
#' @param X matrix of data frame of covariates
#' @param y vector of binary responses
#' @param alpha the alpha as in glmnet. Defauly is 1 (= Lasso). If NULL, no penalty is applied.
#' 
#' @export
baseline.risk <- function(X, y, alpha = 1){
  
  if(is.null(alpha)){
    
    # fit the model
    model.obj <- stats::glm(y~., 
                            family = stats::binomial(link = "logit"), 
                            data = data.frame(y, X))
    
    # get liner predictor etc.
    lp        <- as.numeric(cbind(1,X) %*% model.obj$coefficients)
    kept.vars <- 1:ncol(X)
    
    coefs            <- matrix(model.obj$coefficients)
    rownames(coefs)  <- c("(Intercept)", colnames(X)[kept.vars])
    colnames(coefs)  <- "Estimate"
    lambda           <- NULL
    
  } else{
    
    X             <- as.matrix(X)
    if(is.null(colnames(X))) colnames(X) <- paste0("X", 1:ncol(X))
    model.obj    <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = alpha)
    lambda        <- model.obj$lambda.min # minimizing lambda (needs to be fixed in second stage)
    coefs.obj     <- glmnet::coef.glmnet(model.obj, s = "lambda.min")
    kept.vars     <- coefs.obj@i
    if(0 %in% kept.vars) kept.vars <- kept.vars[-which(kept.vars == 0)] # intercept is not a variable
    coefs            <- as.matrix(coefs.obj@x)
    rownames(coefs)  <- c("(Intercept)", colnames(X)[kept.vars])
    colnames(coefs)  <- "Estimate"
    X.retained       <- cbind(intercept = 1, X[,kept.vars])
    lp               <- as.numeric(X.retained %*% coefs) # linear predictor
    
  } # IF
  
  
  
  return(list(
    model.obj = model.obj,
    lambda.min = lambda,
    linear.predictor = lp,
    response = stats::plogis(lp),
    coefficients = coefs,
    retained.variables = colnames(X)[kept.vars]
  ))
} # FUN


#' perform risk modeling: second stage 
#' \code{y = g(beta_0 +beta_1 * w + beta_2 * z + beta_2 w * z)}
#' 
#' @param z the `z` in the 2nd stage, which is used as product with w and as a single regressor. 
#' @param y a binary response vector
#' @param w a binary treatment assignment vector
#' @param constant.treatment.effect If TRUE, the interaction z*w is not used as regressor. Default is FALSE.
#' @param intercept logical. Shall an intercept be included? Default is `FALSE`
#' 
#' @export 
risk.model.stage2 <- function(z, y, w, 
                              constant.treatment.effect = FALSE,
                              intercept = TRUE){
  
  # prepare flipped W
  w.rev                  <- ifelse(w == 1, 0, 1)
  
  if(constant.treatment.effect){
    
    # prepare X for second stage...
    X.stage2 <- cbind(w = w, z = z)

    # ... and its counterpart with flipped w
    X.stage2.rev           <- cbind(w = w.rev, z = z)

    # prepare formula
    f <- ifelse(intercept, "y ~ w + z", "y ~ 0 + w + z")
  } else{
    
    # prepare X for second stage...
    X.stage2 <- cbind(w = w, z = z, w.z = w * z)
    
    # ... and its counterpart with flipped w
    X.stage2.rev <- cbind(w = w.rev, z = z, w.z = w.rev * z)
    
    # prepare formula
    f <- ifelse(intercept, "y ~ w + z + w.z", "y ~ 0 + w + z + w.z")
    
  } # IF
  
  
  # fit the model
  mod.stage2 <- stats::glm(formula = stats::formula(f), 
                           family = stats::binomial(link = "logit"), 
                           data = as.data.frame(X.stage2))
  
  # get the responses with the regular w ( = F_logistic(x'beta + z))
  risk.regular.w <- as.numeric(stats::predict.glm(mod.stage2, type = "response"))
  
  # get responses with flipped W
  if(intercept){
    X.temp <- cbind(1, X.stage2.rev)  
  } else{
    X.temp <- X.stage2.rev
  } # IF
  
  risk.flipped.w <- as.numeric(stats::plogis(X.temp %*% mod.stage2$coefficients))
  
  # return
  return(list(mod.stage2 = mod.stage2,
              coefficients = summary(mod.stage2)$coefficients,
              risk.regular.w = risk.regular.w,
              risk.flipped.w = risk.flipped.w,
              z = z))
} # FUN


#' performs risk modeling by penalized logistic regression
#' 
#' @param X design matrix or data frame 
#' @param y a binary response vector
#' @param w a binary treatment assignment vector
#' @param alpha the alpha as in glmnet. Default is 1 (= Lasso)
#' @param lifeyears vector of life years. Default is NULL.
#' @param prediction.timeframe vector of the prediction time frame. Default is NULL.
#' @param z the `z` in the 2nd stage, which is used as product with w and as a single regressor. If NULL, thelinear predictor of stage 1 will be used.
#' @param constant.treatment.effect If TRUE, the interaction z*w is not used as regressor in stage 2. Default is FALSE.
#' @param intercept.stage.2 logical. Shall an intercept in stage 2 be included? Default is `TRUE`
#' 
#' @export
risk.modeling <- function(X, y, w, alpha = 1, 
                          z = NULL,
                          lifeyears = NULL,
                          prediction.timeframe = NULL,
                          intercept.stage.2 = TRUE,
                          constant.treatment.effect = FALSE){
  
  # truncate y if necessary
  y.orig    <- y
  
  if(!is.null(lifeyears) & !is.null(prediction.timeframe)){
    lifeyears <- ifelse(lifeyears <= prediction.timeframe, lifeyears, prediction.timeframe) 
    y         <- ifelse(lifeyears <= prediction.timeframe, y, 0)
  } # IF
  
  # make X a matrix
  X <- as.matrix(X)
  
  ## stage 1
  if(is.null(z)){
    stage1 <- baseline.risk(X = X, y = y, alpha = alpha)
    z <- stage1$linear.predictor
  } else{
    stage1 <- NULL
  } # IF
  
  
  ## stage 2
  stage2 <- risk.model.stage2(z = z, y = y, w = w,
                              constant.treatment.effect = constant.treatment.effect, 
                              intercept = intercept.stage.2)
  
  # absolute predicted benefit
  pred.ben.abs.raw <- stage2$risk.regular.w - stage2$risk.flipped.w
  pred.ben.abs     <- ifelse(w == 1, pred.ben.abs.raw, -pred.ben.abs.raw) 
  
  # relative predicted benefit
  pred.ben.rel.raw <- stage2$risk.regular.w / stage2$risk.flipped.w
  pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)
  
  # return
  return(list(
    inputs = list(X = X, w = w, y = y.orig, 
                  lifeyears = lifeyears, 
                  prediction.timeframe = prediction.timeframe, 
                  y.prediction.timeframe = y),
    models = list(model.stage1 = stage1$model.obj,
                  model.stage2 = stage2$mod.stage2,
                  coefficients.stage1 = stage1$coefficients,
                  coefficients.stage2 = stage2$coefficients),
    average.treatment.effect = mean(pred.ben.abs),
    risk = list(risk.regular.w = stage2$risk.regular.w,
                risk.flipped.w = stage2$risk.flipped.w,
                risk.baseline = stage1$response),
    benefits = list(predicted.absolute.benefit = pred.ben.abs,
                    predicted.relative.benefit = pred.ben.rel,
                    predicted.absolute.benefit.raw = pred.ben.abs.raw,
                    predicted.relative.benefit.raw = pred.ben.rel.raw),
    C.statistics = list(c.index.outcome.stage1 = C.index.outcome(y = y, risk.prediction = stage1$response),
                        c.index.outcome.stage2 = C.index.outcome(y = y, risk.prediction = stage2$risk.regular.w),
                        c.index.benefit = C.index.benefit(y = y, w = w, predicted.benefit = pred.ben.abs)),
    z = stage2$z
  ))
} # FUN
