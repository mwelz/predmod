source(paste0(getwd(),  "/funs/c-statistics/c-statistics.R")) # for C index calculations

#' baseline risk prediction via penalized Cox regression (treatment assignment is not used)
#' 
#' @param X matrix of data frame of covariates
#' @param y vector of binary responses
#' @param alpha the alpha as in glmnet. Default is 1 (= Lasso)
#' @param lifeyears vector of life years.
#' @param prediction.timeframe vector of the prediction time frame.
#' 
#' @export
baseline.risk.cox <- function(X, y, 
                              lifeyears, prediction.timeframe, 
                              alpha = 1){
  
  X <- as.matrix(X)
  
  # get survival matrix; to be passed as response of glmnet
  survival.matrix           <- survival::Surv(lifeyears, y)
  colnames(survival.matrix) <- c("time", "status")
  
  # regularized Cox regression for survival
  glmnet.obj <- glmnet::cv.glmnet(x = X, y = survival.matrix, 
                                  family = "cox", type.measure = "C", alpha = alpha)
  lambda     <- glmnet.obj$lambda.min # best lambda (needs to be fixed in second stage)
  coefs.obj  <- glmnet::coef.glmnet(glmnet.obj, s = "lambda.min")
  kept.vars  <- coefs.obj@i + 1 # Cox has no constant, so account for zero-indexing
  coefs      <- coefs.obj@x
  lp         <- as.numeric(predict(glmnet.obj$glmnet.fit, s = lambda,
                           newx = X[,kept.vars],
                           type = "link")) # linear predictors
  
  # base hazard
  basehaz <- hdnom::glmnet_basesurv(time = lifeyears, 
                                    event = y, lp = lp ,
                                    times.eval = prediction.timeframe, 
                                    centered = FALSE)
  
  return(list(
    glmnet.obj = glmnet.obj,
    lambda.min = lambda,
    linear.predictor = lp,
    response = plogis(lp),
    coefficients = coefs,
    retained.variables = colnames(X)[kept.vars],
    timeframe = basehaz$times,
    cumulative.base.hazard = basehaz$cumulative_base_hazard,
    cumulative.survival = basehaz$base_surv
  ))
} # FUN


#' perform risk modeling: second stage (penalized Cox regression):
#' y = g(\beta_0 +\beta_1 * w + \beta_2 w * z + z + \varepsilon)
#' 
#' @param linear.predictor a vector of linear predictions
#' @param y a binary response vector
#' @param w a binary treatment assignment vector
#' @param z the `z` in the model, which is used as product with w and as an offset. Default is `linear.predictor`.
#' @param constant.treatment.effect If TRUE, the interaction z*w is not used as regressor. Default is FALSE.
#' @param lifeyears vector of life years.
#' @param prediction.timeframe vector of the prediction time frame.
#' 
#' @export 
cox.risk.model.stage2  <- function(linear.predictor, y, w, z, 
                                   constant.treatment.effect = FALSE,
                                   lifeyears, prediction.timeframe){
  
  # check input
  if(any(is.character(z))){
    
    if(z == "linear.predictor"){
      zz <- linear.predictor
    } else{
      stop("z needs to be numeric or equal to 'linear predictor'!")
    } # IF
  } else{
    
    zz <- z
    
  } # IF
  
  # prepare flipped W
  w.rev                  <- ifelse(w == 1, 0, 1)
  
  if(constant.treatment.effect){
    
    # prepare X for second stage...
    X.stage2 <- as.matrix(w)
    colnames(X.stage2) <- "w"
    
    # ... and its counterpart with flipped w
    X.stage2.rev           <- as.matrix(w)
    colnames(X.stage2.rev) <- "w"
    
  } else{
    
    # prepare X for second stage...
    X.stage2 <- cbind(w = w, w.z = w * zz)
    
    # ... and its counterpart with flipped w
    X.stage2.rev <- cbind(w = w.rev, w.z = w.rev * zz)
    
  } # IF
  
  # fit the model (no offset due to lack of constant)
  mod.stage2 <- survival::coxph(y ~., data = data.frame(X.stage2, y = survival::Surv(lifeyears, y)))
  
  # get the linear predictor of stage 2 (TODO: does this really account for the centering in the object output?)
  lp.stage2     <- as.numeric(X.stage2 %*% mod.stage2$coefficients)
  lp.stage2.rev <- as.numeric(X.stage2.rev %*% mod.stage2$coefficients)
  
  # base hazard
  stage2.basehaz <- hdnom::glmnet_basesurv(time = lifeyears, event = y, 
                                           lp = lp.stage2, times.eval = prediction.timeframe,
                                           centered = FALSE)
  
  # get the responses with the regular w
  risk.regular.w <- 1 - exp(-stage2.basehaz$cumulative_base_hazard)^exp(lp.stage2)
  
  # get the responses with the flipped w
  risk.flipped.w <- 1 - exp(-stage2.basehaz$cumulative_base_hazard)^exp(lp.stage2.rev)
  
  # return
  return(list(mod.stage2 = mod.stage2,
              coefficients = summary(mod.stage2)$coefficients,
              risk.regular.w = risk.regular.w,
              risk.flipped.w = risk.flipped.w,
              z = zz,
              timeframe = stage2.basehaz$times,
              cumulative.base.hazard = stage2.basehaz$cumulative_base_hazard,
              cumulative.survival = stage2.basehaz$base_surv))
} # FUN



#' performs Cox risk modeling by penalized logistic regression
#' 
#' @param X design matrix or data frame 
#' @param y a binary response vector
#' @param w a binary treatment assignment vector
#' @param alpha the alpha as in glmnet. Default is 1 (= Lasso)
#' @param constant.treatment.effect If TRUE, the interaction z*w is not used as regressor in stage 2. Default is FALSE.
#' @param lifeyears vector of life years. Default is NULL.
#' @param prediction.timeframe vector of the prediction time frame. Default is NULL.
#' @param z the `z` in the 2nd stage, which is used as product with w and as an offset. Default is `linear.predictor`.
#' 
#' @export
cox.risk.modeling <- function(X, w, y, alpha = 1, 
                              constant.treatment.effect = FALSE,
                              lifeyears, prediction.timeframe, z = "linear.predictor"){
  
  # truncate y if necessary
  y.orig    <- y
  lifeyears <- ifelse(lifeyears <= prediction.timeframe, lifeyears, prediction.timeframe) 
  y         <- ifelse(lifeyears <= prediction.timeframe, y, 0)
  
  # make X a matrix
  X <- as.matrix(X)
  
  # stage 1
  stage1 <- baseline.risk.cox(X = X, y = y, lifeyears = lifeyears,
                              prediction.timeframe = prediction.timeframe, alpha = alpha)
  
  # stage 2
  stage2 <- cox.risk.model.stage2(linear.predictor = stage1$linear.predictor, 
                                  y = y, w = w, z = z, 
                                  constant.treatment.effect = constant.treatment.effect, 
                                  lifeyears = lifeyears,
                                  prediction.timeframe = prediction.timeframe)
  
  # absolute predicted benefit
  pred.ben.abs.raw <- stage2$risk.regular.w - stage2$risk.flipped.w
  pred.ben.abs     <- ifelse(w == 1, pred.ben.abs.raw, -pred.ben.abs.raw) 
  
  # relative predicted benefit
  pred.ben.rel.raw <- stage2$risk.regular.w / stage2$risk.flipped.w
  pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)
  
  # C-stat belonging to lambda.min
  lambda.min.index <- which(stage1$glmnet.obj[["lambda"]] %in% stage1$glmnet.obj$lambda.min)
  c.index.outcome.stage1 <- stage1$glmnet.obj$cvm[lambda.min.index]
  
  # return
  return(list(
    inputs = list(X = X, w = w, y = y.orig, 
                  lifeyears = lifeyears, 
                  prediction.timeframe = prediction.timeframe, 
                  y.prediction.timeframe = y),
    models = list(model.stage1 = stage1$glmnet.obj,
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
    C.statistics = list(c.index.outcome.stage1 = c.index.outcome.stage1,
                        c.index.outcome.stage2 = list(estimate = unname(stage2$mod.stage2$concordance["concordance"]),
                                                      stderr   = unname(stage2$mod.stage2$concordance["std"])),
                        c.index.benefit = C.index.benefit(y = y, w = w, predicted.benefit = pred.ben.abs)),
    linear.predictor = stage1$linear.predictor,
    z = stage2$z
  ))
} # FUN
