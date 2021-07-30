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
#' @param lifeyears vector of life years.
#' @param prediction.timeframe vector of the prediction time frame.
#' @param lambda the lambda for the penalty term
#' 
#' @export 
cox.risk.model.stage2  <- function(linear.predictor, y, w, z, 
                                   lifeyears, prediction.timeframe, lambda){
  
  # check input for offset.linear.predictor
  if(any(is.character(z))){
    
    if(z == "linear.predictor"){
      z <- linear.predictor
    } else{
      stop("z needs to be numeric or equal to 'linear predictor'!")
    } # IF
  } # IF
  
  # prepare X for second stage
  X.stage2 <- cbind(w = w, w.z = w * z)
  
  # get survival matrix; to be passed as response of glmnet
  survival.matrix           <-survival::Surv(lifeyears, y)
  colnames(survival.matrix) <- c("time", "status")
  
  # Cox regression for survival
  mod.stage2 <- glmnet::glmnet(x = X.stage2, y = survival.matrix,
                               family = "cox", type.measure = "C",
                               lambda = lambda,
                               offset = z,
                               alpha = alpha) # TODO: ask Kevin, no offset previously!
  
  # get the estimated coefficients
  coefs.obj     <- glmnet::coef.glmnet(mod.stage2)
  coefs         <- coefs.obj@x
  kept.vars     <- coefs.obj@i + 1 # Cox has no constant, so account for zero-indexing
  
  # prepare design matrix with retained variables
  X.retained  <- X.stage2[,kept.vars]
  
  # get the linear predictor of stage 2 (including the offset)
  lp.stage2 <- as.numeric(X.retained %*% coefs) + as.numeric(z)
  
  # base hazard
  stage2.basehaz <- hdnom::glmnet_basesurv(time = lifeyears, event = y, 
                                           lp = lp.stage2, times.eval = prediction.timeframe,
                                           centered = FALSE)
  
  # get the responses with the regular w
  risk.regular.w <- 1 - exp(-stage2.basehaz$cumulative_base_hazard)^exp(lp.stage2)
  
  # get the responses with flipped w
  w.rev           <- ifelse(w == 1, 0, 1)
  X.stage2.rev    <- cbind(w = w.rev, w.z = w.rev * z)
  
  # prepare design matrix with flipped w
  X.retained.rev  <- X.stage2.rev[,kept.vars]
  
  # get the linear predictor of stage 2 with flipped w (inc. offset)
  lp.stage2.rev <- as.numeric(X.retained.rev %*% coefs) + as.numeric(z)
  
  # calculate risk with flipped w
  risk.flipped.w <-  1 - exp(-stage2.basehaz$cumulative_base_hazard)^exp(lp.stage2.rev)
  
  # return
  return(list(mod.stage2 = mod.stage2,
              risk.regular.w = risk.regular.w,
              risk.flipped.w = risk.flipped.w,
              z = z,
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
#' @param lifeyears vector of life years. Default is NULL.
#' @param prediction.timeframe vector of the prediction time frame. Default is NULL.
#' @param z the `z` in the 2nd stage, which is used as product with w and as an offset. Default is `linear.predictor`.
#' 
#' @export
cox.risk.modeling <- function(X, w, y, alpha = 1, 
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
                                  lifeyears = lifeyears, prediction.timeframe = prediction.timeframe, 
                                  lambda = stage1$lambda.min)
  
  # absolute predicted benefit
  pred.ben.abs.raw <- stage2$risk.regular.w - stage2$risk.flipped.w
  pred.ben.abs     <- ifelse(w == 1, pred.ben.abs.raw, -pred.ben.abs.raw) 
  
  # relative predicted benefit
  pred.ben.rel.raw <- stage2$risk.regular.w / stage2$risk.flipped.w
  pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)
  
  # coefficients
  coefs.stage1 <- as.matrix(glmnet::coef.glmnet(stage1$glmnet.obj, s = "lambda.min"))
  coefs.stage2 <- as.matrix(glmnet::coef.glmnet(stage2$mod.stage2))
  colnames(coefs.stage2) <- colnames(coefs.stage1) <- "Estimated Coefficient"
  
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
                  coefficients.stage1 = coefs.stage1,
                  coefficients.stage2 = coefs.stage2),
    average.treatment.effect = mean(pred.ben.abs),
    risk = list(risk.regular.w = stage2$risk.regular.w,
                risk.flipped.w = stage2$risk.flipped.w,
                risk.baseline = stage1$response),
    benefits = list(predicted.absolute.benefit = pred.ben.abs,
                    predicted.relative.benefit = pred.ben.rel,
                    predicted.absolute.benefit.raw = pred.ben.abs.raw,
                    predicted.relative.benefit.raw = pred.ben.rel.raw),
    C.statistics = list(c.index.outcome.stage1 = c.index.outcome.stage1,
                        c.index.outcome.stage2 = C.index.outcome(y = y, risk.prediction = stage2$risk.regular.w),
                        c.index.benefit = C.index.benefit(y = y, w = w, predicted.benefit = pred.ben.abs)),
    linear.predictor = stage1$linear.predictor,
    z = stage2$z
  ))
} # FUN
