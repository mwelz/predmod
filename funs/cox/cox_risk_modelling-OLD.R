cox.risk.model.stage1 <- function(X, y,lifeyears, predictiontimeframe, alpha = 1){
  
  # stage 1 doesn't use W!
  # in 

  
  glmsurvival.coxstage1.obj <-survival::Surv(lifeyears, y)
  colnames(glmsurvival.coxstage1.obj) <- c("time", "status")
  cox.models.stage1 <- glmnet::cv.glmnet(x, glmsurvival.coxstage1.obj, family = "cox", type.measure = "C", alpha = 1)
  lambda        <- cox.models.stage1$lambda.min # best lambda (needs to be fixed in second stage)
  coefs.obj     <- glmnet::coef.glmnet(cox.models.stage1, s = "lambda.min")
  kept.vars     <- (coefs.obj@i)+1 #cox has no constant; ignore centering value
 # if(0 %in% kept.vars) kept.vars <- kept.vars[-which(kept.vars == 0)]
  coefs         <- coefs.obj@x
  X.star        <- cbind(x[,kept.vars])
  lp            <-predict(cox.models.stage1$glmnet.fit, s=lambda, X.star,type="link") #linear predictors
  #lp            <- predict(cox.models.stage1, newx = X, s = "lambda.min") #as.numeric(X.star %*% coefs)
  stage1.basehaz <- hdnom::glmnet_basesurv(lifeyears, y, lp , times.eval = predictiontimeframe, centered = FALSE)
  #lp2            <- as.numeric(X.star %*% coefs) check for equality to predict(); provides the sme linear predictors
 
  return(list(
   cox.models.stage1.cv = cox.models.stage1,
    lambda.min = lambda,
    lp = lp,
   stage1.timeframe = stage1.basehaz$times,
   stage1.cum.basehaz = stage1.basehaz$cumulative_base_hazard,
   stage1.cum.surv = stage1.basehaz$base_surv
  ))
}


cox.risk.model.stage2  <- function(lp, y, lifeyears, w, predictiontimeframe, lambda, offset.lp = TRUE ){
  
  #X.stage2 <- cbind(w = w, wlp = w * lp, lp)
  X.stage2 <- cbind(w = w, lp= lp)
  if(offset.lp){
    oset <- lp
  } else{
    oset <- rep(0, length(lp))
  }
  glmsurvival.coxstage2.obj <-survival::Surv(lifeyears, y)
  colnames(glmsurvival.coxstage2.obj) <- c("time", "status")
  # use same lambda as model 1!
  mod.stage2 <- glmnet::glmnet(X.stage2, glmsurvival.coxstage2.obj, family = "cox", type.measure = "C", lambda = lambda)
  coefs          <- mod.stage2$beta
 # lp.stage2     <-predict(mod.stage2, s=lambda, X.stage2,type="link") #linear predictors
  lp.stage2      <- as.numeric(X.stage2  %*% coefs)
  stage2.basehaz <- hdnom::glmnet_basesurv(lifeyears, y, lp.stage2  , predictiontimeframe, centered = FALSE)
  
  
  ### risk prediction:
  ## Note: This might be a bug in glmnet, but a nonzero offset causes the predicted 'probabilities' to 
  # not be bounded by [0,1] anymore. Hence we need to apply to logistic transformation on them one more time
  # (we do so in the 'return' line)
  

  # 1) ordinary w
  probs.regular <-  1-exp(-stage2.basehaz$cumulative_base_hazard)^exp(lp.stage2)
  # 2) w flipped
  w.rev           <- ifelse(w == 1, 0, 1)
  X.stage2.rev    <- cbind(w = w.rev, wlp = w.rev * lp)
  lp.stage2.rev   <- as.numeric(X.stage2.rev  %*% coefs)
  probs.flipped.w <-  1-exp(-stage2.basehaz$cumulative_base_hazard)^exp(lp.stage2.rev)
  
  # return
  return(list(mod.stage2 = mod.stage2,
              risk.pred.regular = (as.numeric(probs.regular)),
              risk.pred.flipped.w = (as.numeric(probs.flipped.w)),
              stage2.timeframe = stage2.basehaz$times,
              stage2.cum.basehaz = stage2.basehaz$cumulative_base_hazard,
              stage2.cum.surv = stage2.basehaz$base_surv
  ))
}


cox.risk.modeling <- function(X, w, y, lifeyears, predictiontimeframe, alpha, offset.lp = TRUE){
  ## stage 1
  lifeyears <- ifelse(lifeyears <=predictiontimeframe, lifeyears, predictiontimeframe)
  y<- ifelse(lifeyears <=predictiontimeframe, y, 0)
  stage1 <- cox.risk.model.stage1(X = X, y = y, lifeyears = lifeyears, predictiontimeframe = predictiontimeframe, alpha = alpha)

  ## stage 2
  stage2 <- cox.risk.model.stage2(lp = stage1$lp,
                              y = y, 
                              lifeyears = lifeyears,
                              w = w,
                              predictiontimeframe = predictiontimeframe,
                              lambda = stage1$lambda.min,
                              offset.lp = offset.lp)
  
  
  # absolute predicted benefit
  pred.ben.abs.raw <- stage2$risk.pred.regular - stage2$risk.pred.flipped.w
  pred.ben.abs     <- ifelse(w == 1, pred.ben.abs.raw, -pred.ben.abs.raw) 
  
  # relative predicted benefit
  pred.ben.rel.raw <- stage2$risk.pred.regular / stage2$risk.pred.flipped.w
  pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)

  # coefficients
  coefs.stage1 <-  as.matrix(coef(stage1$cox.models.stage1.cv, s = "lambda.min"))
  coefs.stage2 <- as.matrix(glmnet::coef.glmnet(stage2$mod.stage2))
  colnames(coefs.stage2) <- colnames(coefs.stage1) <- "Estimated Coefficient"
  
  #Match cases based on observed benefit
  Treatment.formula<- w~pred.ben.abs
  matched <- MatchIt::matchit(Treatment.formula)
  treated <- as.numeric(rownames(matched$match.matrix))
  control <- as.numeric(matched$match.matrix[,1])
  
  #Remove unpaired observations
  no.pairing <- which(is.na(matched$match.matrix))
  treated <- treated[-no.pairing]
  control <- control[-no.pairing]
  
  obs.ben <- y[control]-y[treated]
  pred.ben.abs.paired = (pred.ben.abs[control]+pred.ben.abs[treated])/2 #Pairs are matched by predicted benefit; average over the pair
  
  # calculate C for benefit by using predicted risk (with regular w)
  c.index.benefit = unname(Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)[1])
  c.index.stage1.youtcome <- unname(Hmisc::rcorr.cens(stage1$lp, y)[1])
  pred.ben.abs.paired = pred.ben.abs[control]-pred.ben.abs[treated]
  c.index.benefit = unname(Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)[1])
  
  # calculate C for benefit by using predicted risk (with regular w)
  #C-stat belonging to lambda.min
  lambda.min.index = which(stage1$cox.models.stage1.cv[["lambda"]] %in% stage1$cox.models.stage1.cv$lambda.min)
  c.index.youtcome = stage1$cox.models.stage1.cv[["cvm"]][lambda.min.index] 
  
  return(list(
    inputs = list(X = X, w = w, y = y),
    mod.stage1 = stage1$models.stage1.cv,
    mod.stage2 = stage2$mod.stage2,
    coefficients.stage1 = coefs.stage1,
    coefficients.stage2 = coefs.stage2,
    linear.predictor = stage1$lp,
    risk.baseline =  1-exp(-stage1$stage1.cum.basehaz)^exp(stage1$lp),
    risk.regular.w = stage2$risk.pred.regular,
    risk.flipped.w = stage2$risk.pred.flipped.w,
    predicted.absolute.benefit = pred.ben.abs,
    predicted.absolute.benefit.raw = pred.ben.abs.raw,
    predicted.relative.benefit = pred.ben.rel,
    predicted.relative.benefit.raw = pred.ben.rel.raw,
    ate.hat = mean(pred.ben.abs),
    c.index.benefit = c.index.benefit,
    c.index.youtcome = c.index.youtcome
  ))
}
