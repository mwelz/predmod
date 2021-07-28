# load imputation accounters
source(paste0(getwd(), "/funs/imputation/imputation.R"))

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
  names(coefs)  <- c("(Intercept)", colnames(X)[kept.vars])
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


#' perform risk modeling: second stage (penalized logistic regression):
#' y = g(\beta_0 +\beta_1 * w + \beta_2 w * z + z + \varepsilon)
#' 
#' @param linear.predictor a vector of linear predictions
#' @param y a binary response vector
#' @param w a binary treatment assignment vector
#' @param z the `z` in the model, which is used as product with w and as an offset. Default is `linear.predictor`.
#' @param lambda the lambda for the penalty term
#' @param intercept logical. Shall an intercept be included? Default is `FALSE`
#' 
#' @export 
risk.model.stage2 <- function(linear.predictor, y, w, z, lambda, 
                              intercept = FALSE){
  
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
  
  # stage 2 modeling
  mod.stage2 <- glmnet::glmnet(X.stage2, y, 
                               family = "binomial",
                               lambda = lambda,
                               intercept = intercept,
                               offset = z) 
  
  # get the estimated coefficients
  coefs.obj     <- glmnet::coef.glmnet(mod.stage2)
  kept.vars     <- coefs.obj@i
  if(0 %in% kept.vars) kept.vars <- kept.vars[-which(kept.vars == 0)] # intercept is not a variable
  coefs         <- coefs.obj@x
  
  # prepare design matrix with retained variables
  if(intercept){
    X.retained  <- cbind(intercept = 1, X.stage2[,kept.vars])
  } else{
    X.retained  <- X.stage2[,kept.vars]
  } # IF
  
  # get the responses with the regular w
  risk.regular.w <- plogis(as.numeric(X.retained %*% coefs) + as.numeric(z))
  
  # get the responses with flipped w
  w.rev           <- ifelse(w == 1, 0, 1)
  X.stage2.rev    <- cbind(w = w.rev, w.z = w.rev * z)
  
  # prepare design matrix with flipped w
  if(intercept){
    X.retained.rev  <- cbind(intercept = 1, X.stage2.rev[,kept.vars])
  } else{
    X.retained.rev  <- X.stage2.rev[,kept.vars]
  } # IF
  
  # calculate risk with flipped w
  risk.flipped.w  <- plogis(as.numeric(X.retained.rev %*% coefs) + as.numeric(z))
  
  # return
  return(list(mod.stage2 = mod.stage2,
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
#' @param z the `z` in the 2nd stage, which is used as product with w and as an offset. Default is `linear.predictor`.
#' @param lambda the lambda for the penalty term
#' @param intercept.stage.2 logical. Shall an intercept in stage 2 be included? Default is `FALSE`
#' 
#' @export
risk.modeling <- function(X, y, w, alpha = 1, 
                          lifeyears = NULL,
                          prediction.timeframe = NULL,
                          intercept.stage.2 = FALSE,
                          z = "linear.predictor"){
  
  # truncate y if necessary
  y.orig    <- y
  
  if(!is.null(lifeyears) & !is.null(prediction.timeframe)){
    lifeyears <- ifelse(lifeyears <= prediction.timeframe, lifeyears, prediction.timeframe) 
    y         <- ifelse(lifeyears <= prediction.timeframe, y, 0)
  } # IF
  
  # make X a matrix
  X <- as.matrix(X)
  
  ## stage 1
  stage1 <- baseline.risk(X = X, y = y, alpha = alpha)
  
  ## stage 2
  stage2 <- risk.model.stage2(linear.predictor = stage1$linear.predictor,
                              y = y, w = w,
                              lambda = stage1$lambda.min, 
                              intercept = intercept.stage.2,
                              z = z)
  
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
  matched <- MatchIt::matchit(w ~ pb, data = data.frame(w=w, pb=pred.ben.abs))
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
                        c.index.outcome.stage2 = c.index.outcome.stage2,
                        c.index.benefit = c.index.benefit),
    linear.predictor = stage1$linear.predictor,
    z = stage2$z
  ))
} # FUN


#' TODO: write documentation
#'
#'
risk.modeling_imputation.accounter <- function(predictive.model.imputed){
  
  # initialize
  pred.model.imp.adj <- list()
  m <- length(predictive.model.imputed)
  
  # stage 1 coefficients
  pred.model.imp.adj$models$coefficients.stage1 <- imputation.accounter_location(
    lapply(1:m, function(i){
      
      # get names of all variables (pre-selection)
      nam.all.variables <<- c("(Intercept)", colnames(predictive.model.imputed[[i]]$inputs$X))
      
      # initialize long array with zeros for unselected variables
      selected.variables.long <<- rep(0.0, length(nam.all.variables))
      names(selected.variables.long) <<- nam.all.variables
      
      # assign values to long array
      selected.variables.short <<- predictive.model.imputed[[i]]$models$coefficients.stage1
      selected.variables.long[names(selected.variables.short)] <<- selected.variables.short
      selected.variables.long
    }))
  
  # stage 2 coefficients
  pred.model.imp.adj$models$coefficients.stage2 <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$models$coefficients.stage2))
  
  # ATE
  pred.model.imp.adj$average.treatment.effect <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$average.treatment.effect))
  
  # risk regular w
  pred.model.imp.adj$risk$risk.regular.w <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$risk$risk.regular.w))
  
  # risk flipped w
  pred.model.imp.adj$risk$risk.flipped.w <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$risk$risk.flipped.w))
  
  # risk baseline
  pred.model.imp.adj$risk$risk.baseline <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$risk$risk.baseline))
  
  # predicted absolute benefit
  pred.model.imp.adj$benefits$predicted.absolute.benefit <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$benefits$predicted.absolute.benefit))
  
  # predicted relative benefit
  pred.model.imp.adj$benefits$predicted.relative.benefit <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$benefits$predicted.relative.benefit))
  
  # predicted absolute benefit raw
  pred.model.imp.adj$benefits$predicted.absolute.benefit.raw <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$benefits$predicted.absolute.benefit.raw))
  
  # predicted relative benefit raw
  pred.model.imp.adj$benefits$predicted.relative.benefit.raw <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$benefits$predicted.relative.benefit.raw))
  
  # C stat stage 1
  pred.model.imp.adj$C.statistics$c.index.outcome.stage1 <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$C.statistics$c.index.outcome.stage1))
  
  # C stat stage 2
  pred.model.imp.adj$C.statistics$c.index.outcome.stage2 <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$C.statistics$c.index.outcome.stage2))
  
  # C index benefit
  pred.model.imp.adj$C.statistics$c.index.benefit <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$C.statistics$c.index.benefit))
  
  # LP
  pred.model.imp.adj$linear.predictor <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$linear.predictor))
  
  # z
  pred.model.imp.adj$z <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$z))
  
  # return
  return(pred.model.imp.adj)
  
} # FUN