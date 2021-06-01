
risk.model.poiss.stage1 <- function(X, y,ly, alpha = 1){
  
  # stage 1 doesn't use W!
  models.stage1 <- glmnet::cv.glmnet(X, y,offset = log(ly), family = poisson(link = "log"), alpha = alpha)
  lambda        <- models.stage1$lambda.min # best lambda (needs to be fixed in second stage)
  coefs.obj     <- glmnet::coef.glmnet(models.stage1, s = "lambda.min")
  kept.vars     <- coefs.obj@i
  if(0 %in% kept.vars) kept.vars <- kept.vars[-which(kept.vars == 0)]
  coefs         <- coefs.obj@x
  X.star        <- cbind(intercept = 1, X[,kept.vars])
  lp_raw        <- as.numeric(exp(X.star %*% coefs))
  lp            <- lp_raw*ly
  
  return(list(
    models.stage1.cv = models.stage1,
    lambda.min = lambda,
    lp = lp
  ))
}


risk.model.poiss.stage2 <- function(lp , y, w, lambda, offset.lp  = TRUE ){
  
  X.stage2 <- cbind(w = w, wlp = w * lp )
  if(offset.lp){
    oset <- lp
  } else{
    oset <- rep(0, length(lp))
  }
  
  
  # use same lambda as model 1!
  mod.stage2 <- glmnet::glmnet(X.stage2, y, 
                               family = poisson(link = "log"),
                               lambda = lambda,
                               intercept = TRUE,
                               offset = oset) 
  
  ### risk prediction:
  ## Note: This might be a bug in glmnet, but a nonzero offset causes the predicted 'probabilities' to 
  # not be bounded by [0,1] anymore. Hence we need to apply to logistic transformation on them one more time
  # (we do so in the 'return' line)
  
  # 1) ordinary w
  probs.regular <- glmnet::predict.glmnet(mod.stage2, 
                                          newx = X.stage2,
                                          type = "response", 
                                          newoffset = oset)
  probs.regular <- exp(probs.regular)
  
  
  # 2) w flipped
  w.rev           <- ifelse(w == 1, 0, 1)
  X.stage2.rev    <- cbind(w = w.rev, wlp = w.rev * lp)
  probs.flipped.w <- glmnet::predict.glmnet(mod.stage2, 
                                            newx = X.stage2.rev,
                                            type = "response", 
                                            newoffset = oset)
  
  probs.flipped.w <- exp(probs.flipped.w)
  
  # return
  return(list(mod.stage2 = mod.stage2,
              risk.pred.regular = (as.numeric(probs.regular)),
              risk.pred.flipped.w = (as.numeric(probs.flipped.w))
  ))
}


risk.modeling.poiss <- function(X, w, y, ly, alpha, offset.lp = TRUE){
  ## stage 1
  stage1 <- risk.model.poiss.stage1(X = X, y = y, ly=ly, alpha = alpha)
  
  ## stage 2
  stage2 <- risk.model.poiss.stage2(lp = stage1$lp,
                                    y = y, 
                                    w = w,
                                    lambda = stage1$lambda.min,
                                    offset.lp = offset.lp)
  
  # absolute predicted benefit
  pred.ben.abs.raw <- stage2$risk.pred.regular - stage2$risk.pred.flipped.w
  pred.ben.abs     <- ifelse(w == 1, -pred.ben.abs.raw, pred.ben.abs.raw)
  
  Estimated_reduction_per_1000_lY =  (sum(pred.ben.abs)/sum(ly))*1000
  
  # relative predicted benefit
  pred.ben.rel.raw <- stage2$risk.pred.regular / stage2$risk.pred.flipped.w
  pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)
  
  # coefficients
  coefs.stage1 <- as.matrix(glmnet::coef.glmnet(stage1$models.stage1.cv, s = "lambda.min"))
  coefs.stage2 <- as.matrix(glmnet::coef.glmnet(stage2$mod.stage2))
  colnames(coefs.stage2) <- colnames(coefs.stage1) <- "Estimated Coefficient"
  
  # calculate C index by using predicted risk (with regular w)
  c.index <- unname(Hmisc::rcorr.cens(x = stage2$risk.pred.regular, S = y)[1])
  
  return(list(
    inputs = list(X = X, w = w, y = y),
    mod.stage1 = stage1$models.stage1.cv,
    mod.stage2 = stage2$mod.stage2,
    coefficients.stage1 = coefs.stage1,
    coefficients.stage2 = coefs.stage2,
    poiss.baseline = stage1$lp,
    poiss.stage2.regular.w = stage2$risk.pred.regular,
    poiss.stage2.flipped.w = stage2$risk.pred.flipped.w,
    predicted.absolute.benefit = pred.ben.abs,
    predicted.absolute.benefit.raw = pred.ben.abs.raw,
    predicted.relative.benefit = pred.ben.rel,
    predicted.relative.benefit.raw = pred.ben.rel.raw,
    ate.hat = mean(pred.ben.abs),
    rel.hat = mean(pred.ben.rel),
    ratereduction.1000ly = Estimated_reduction_per_1000_lY,
    c.index = NA
    #c.index = c.index
  ))
}

