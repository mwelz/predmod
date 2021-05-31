# TODO: incorporate lifeyears and predictiontimeframe as arguments
# TODO: make output nicer everywhere


#' baseline risk prediction via penalized logistic regression (treatment assignment is not used)
#' 
#' @param X matrix of data frame of covariates
#' @param y vector of binary responses
#' @param alpha the alpha as in glmnet. Defauly is 1 (= Lasso)
#' 
#' @export
baseline.risk <- function(X, y, alpha = 1){
  
  X             <- as.matrix(X)
  if(is.null(colnames(X))) colnames(X) <- paste0("X", 1:ncol(X))
  glmnet.obj    <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = alpha)
  lambda        <- glmnet.obj$lambda.min # minimizing lambda (needs to be fixed in second stage)
  coefs.obj     <- glmnet::coef.glmnet(glmnet.obj, s = "lambda.min")
  kept.vars     <- coefs.obj@i
  if(0 %in% kept.vars) kept.vars <- kept.vars[-which(kept.vars == 0)] # intercept is not a variable
  coefs         <- coefs.obj@x
  names(coefs)  <- c("(Intercept)", colnames(X))
  X.retained    <- cbind(intercept = 1, X[,kept.vars])
  lp            <- as.numeric(X.retained %*% coefs) # linear predictor
  
  return(list(
    glmnet.obj = glmnet.obj,
    lambda.min = lambda,
    linear.predictor = lp,
    response = plogis(lp),
    coefficients = coefs,
    retained.variables = colnames(X)[kept.vars]
  ))
} # FUN


#' perform risk modeling: second stage (penalized logistic regression)
#' 
#' @param linear.predictor a vector of linear predictions
#' @param y a binary response vector
#' @param w a binary treatment assignment vector
#' @param lambda the lambda for the penalty term
#' @param offset.linear.predictor handling an offset. If TRUE (default), then _linear.predictor_ is used as offset. If FALSE, no offset is used. If a numeric vector is supplied, this vector is then used as offset.
#' 
#' @export 
risk.model.stage2 <- function(linear.predictor, y, w, lambda, 
                              offset.linear.predictor = TRUE){
  
  # check input for offset.linear.predictor
  if(is.logical(offset.linear.predictor)){
    
    if(offset.linear.predictor){
      offset <- linear.predictor
    } else{
      offset <- rep(0.0, length(linear.predictor))
    } # IF
    
  } else{
    
    offset <- offset.linear.predictor
    
  } # IF
  
  # prepare X for second stage
  X.stage2 <- cbind(w = w, wlp = w * linear.predictor)
  
  # stage 2 modeling
  mod.stage2 <- glmnet::glmnet(X.stage2, y, 
                               family = "binomial",
                               lambda = lambda,
                               intercept = TRUE,
                               offset = offset) 
  
  # get the estimated coefficients
  coefs.obj     <- glmnet::coef.glmnet(mod.stage2)
  kept.vars     <- coefs.obj@i
  if(0 %in% kept.vars) kept.vars <- kept.vars[-which(kept.vars == 0)] # intercept is not a variable
  coefs         <- coefs.obj@x
  X.retained    <- cbind(intercept = 1, X.stage2[,kept.vars])
  
  # get the responses with the regular w
  risk.regular.w <- plogis(as.numeric(X.retained %*% coefs) + as.numeric(offset))
  
  # get the responses with flipped w
  w.rev           <- ifelse(w == 1, 0, 1)
  X.stage2.rev    <- cbind(w = w.rev, wlp = w.rev * linear.predictor)
  X.retained.rev  <- cbind(intercept = 1, X.stage2.rev[,kept.vars])
  risk.flipped.w  <- plogis(as.numeric(X.retained.rev %*% coefs) + as.numeric(offset))
  
  # return
  return(list(mod.stage2 = mod.stage2,
              risk.regular.w = risk.regular.w,
              risk.flipped.w = risk.flipped.w))
} # FUN


#' performs risk modeling by penalized logistic regression
#' 
#' @param X design matrix or data frame 
#' @param y a binary response vector
#' @param w a binary treatment assignment vector
#' @param alpha the alpha as in glmnet. Default is 1 (= Lasso)
#' @param offset.linear.predictor handling an offset. If TRUE (default), then _linear.predictor_ is used as offset. If FALSE, no offset is used. If a numeric vector is supplied, this vector is then used as offset.
#' 
#' @export
risk.modeling <- function(X, y, w, alpha = 1, offset.linear.predictor = TRUE){
  
  # lifeyears <- ifelse(lifeyears <=predictiontimeframe, lifeyears, predictiontimeframe) TODO!
  # y<- ifelse(lifeyears <=predictiontimeframe, y, 0)
  
  # make X a matrix
  X <- as.matrix(X)
  
  ## stage 1
  stage1 <- baseline.risk(X = X, y = y, alpha = alpha)
  
  ## stage 2
  stage2 <- risk.model.stage2(linear.predictor = stage1$linear.predictor,
                              y = y, w = w,
                              lambda = stage1$lambda.min, 
                              offset.linear.predictor = offset.linear.predictor)
  
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
  
  
  # match cases based on observed benefit
  matched <- MatchIt::matchit(w ~ pred.ben.abs)
  match.treated <- as.numeric(rownames(matched$match.matrix))
  match.control <- as.numeric(matched$match.matrix[,1])
  
  # remove unpaired observations
  no.pairing <- which(is.na(match.control))
  if(length(no.pairing) > 0){
    match.treated <- match.treated[-no.pairing]
    match.control <- match.control[-no.pairing]
  }
  
  
  # calculate C for benefit by using predicted risk (with regular w)
  obs.ben             <- y[match.control] - y[match.treated]
  pred.ben.abs.paired <- (pred.ben.abs[match.control] +
                            pred.ben.abs[match.treated]) / 2
  c.index.benefit     <- unname(Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)[1])
  
  # C index
  c.index.outcome.stage1 <- unname(Hmisc::rcorr.cens(stage1$response, y)[1])
  c.index.outcome.stage2 <- unname(Hmisc::rcorr.cens(stage2$risk.regular.w, y)[1])
  
  
  # return
  return(list(
    inputs = list(X = X, w = w, y = y),
    mod.stage1 = stage1$glmnet.obj,
    mod.stage2 = stage2$mod.stage2,
    coefficients.stage1 = coefs.stage1,
    coefficients.stage2 = coefs.stage2,
    linear.predictor = stage1$linear.predictor,
    risk.baseline = stage1$response,
    risk.regular.w = stage2$risk.regular.w,
    risk.flipped.w = stage2$risk.flipped.w,
    predicted.absolute.benefit = pred.ben.abs,
    predicted.absolute.benefit.raw = pred.ben.abs.raw,
    predicted.relative.benefit = pred.ben.rel,
    predicted.relative.benefit.raw = pred.ben.rel.raw,
    ate.hat = mean(pred.ben.abs),
    c.index.benefit = c.index.benefit,
    c.index.outcome.stage1 = c.index.outcome.stage1,
    c.index.outcome.stage2 = c.index.outcome.stage2
  ))
} # FUN

