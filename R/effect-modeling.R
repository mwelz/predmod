#' get the matrix \code{"w * X[, interactions]}"
#' 
#' @param X design matrix
#' @param w vector of binary treatment assignments
#' @param interacted.variables array of strings of the variables in X that shall be interacted with w. If equal to "all", all variables in X are interacted with _w_.
#' 
#' @noRd
get.interaction.terms.matrix <- function(X, w, interacted.variables){
  
  if(is.null(interacted.variables)) return(NULL)
  
  if("w" %in% interacted.variables) stop("w cannot be interacted with itself!")
  
  # add interaction variables
  interaction.terms.matrix <- sapply(which(colnames(X) %in% interacted.variables),
                                     function(j) ifelse(w == 1, X[,j], 0))
  
  colnames(interaction.terms.matrix) <- paste0("w.", interacted.variables)
  
  return(interaction.terms.matrix)
  
} # FUN


# interacted.variables is either NULL or variable names (not indices!)
effect.model.predicted.benefits <- function(final.model, X, kept.vars, interacted.variables){
  
  # predict risk with ordinary W, Pr(Y=1 |W,X)
  response <- unname(stats::predict.glm(final.model, type = "response"))
  
  # flip w
  w <- final.model$data$w
  w.flipped   <- ifelse(w == 1, 0, 1)
  
  # get temporary Z matrix
  Z.temp <- cbind(X, w = w.flipped,
                  get.interaction.terms.matrix(X = X, w = w.flipped, 
                                               interacted.variables = interacted.variables))
  
  # keep Z with flipped W
  Z.flipped <- Z.temp[,kept.vars, drop = FALSE]
  
  # predict risk with flipped W
  response_w.flipped <- unname(stats::predict.glm(final.model,
                                                  newdata = as.data.frame(Z.flipped), 
                                                  type = "response"))
  
  # get absolute predicted benefit
  pred.ben.abs.raw <- response - response_w.flipped
  pred.ben.abs     <- ifelse(w == 1, pred.ben.abs.raw, -pred.ben.abs.raw)
  
  # get relative predicted benefit
  pred.ben.rel.raw <- response / response_w.flipped
  pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)
  
  # return
  return(list(pred.ben.abs.raw = pred.ben.abs.raw, 
              pred.ben.abs = pred.ben.abs,
              pred.ben.rel.raw = pred.ben.rel.raw,
              pred.ben.rel = pred.ben.rel,
              risk.regular.w = response,
              risk.flipped.w = response_w.flipped))
  
} # FUN


#' performs effect modeling by penalized logistic regression
#' 
#' TODO: make argument to pass baseline risk as argument
#' TODO: pass all arguments of baselinerisk() and allow for different alpha
#' 
#' @param X design matrix or data frame 
#' @param y a binary response vector
#' @param w a binary treatment assignment vector
#' @param alpha the alpha as in glmnet. Default is 1 (= Lasso)
#' @param lifeyears vector of life years. Default is NULL.
#' @param prediction.timeframe vector of the prediction time frame. Default is NULL.
#' @param interacted.variables character vector of variables that are to be interacted with \code{w}
#' @param retained character vector of variables that shall not be regularized
#' @param retain_w Logical. Shall \code{w} be retained?
#' 
#' @export
effect.modeling <- function(X, y, w, alpha = 1, 
                            lifeyears = NULL,
                            prediction.timeframe = NULL,
                            interacted.variables = NULL,
                            retained = FALSE,
                            retain_w = TRUE){
  
  # truncate y if necessary
  y.orig    <- y
  
  if(!is.null(lifeyears) & !is.null(prediction.timeframe)){
    lifeyears <- ifelse(lifeyears <= prediction.timeframe, lifeyears, prediction.timeframe) 
    y         <- ifelse(lifeyears <= prediction.timeframe, y, 0)
  } # IF
  
  # make X a matrix
  X <- as.matrix(X)
  
  # get matrix of interaction effects
  Z <- cbind(X, w,
             get.interaction.terms.matrix(X = X, w = w, interacted.variables = interacted.variables))
  
  # get index of W
  idx.w <- which(colnames(Z) %in% "w")
  
  # binary penalty factor. Zero means that corresponding variable isn't shrunk
  penalty.factor <- rep(1, ncol(Z))
  
  if(!is.null(retained)){
    retained.idx   <- which(colnames(Z) %in% retained)
    penalty.factor[retained.idx] <- 0
  } # IF
  
  
  if(retain_w){
    penalty.factor[idx.w] <- 0 
  } # IF
  
  
  # fit penalized regression
  mod <- glmnet::cv.glmnet(x = Z, y = y, 
                           family = "binomial",
                           alpha = alpha, 
                           penalty.factor = penalty.factor)
  
  # extract retained variables
  coefs.obj <- glmnet::coef.glmnet(mod, s = "lambda.min")
  kept.vars <- coefs.obj@i
  
  # intercept isn't a covariate
  if(0 %in% kept.vars){ 
    kept.vars <- kept.vars[-which(kept.vars == 0)]
  } # IF
  
  # fit final model
  Z.retained  <- Z[,kept.vars, drop = FALSE]
  final.model <- stats::glm(y~., data = data.frame(y = y, Z.retained),
                            family =  stats::binomial(link = "logit"))
  
  # calculate benefits
  ben.obj <- effect.model.predicted.benefits(final.model, X, kept.vars, interacted.variables)
  
  # extract benefits
  pred.ben.abs.raw <- ben.obj$pred.ben.abs.raw
  pred.ben.abs <- ben.obj$pred.ben.abs
  pred.ben.rel <- ben.obj$pred.ben.rel
  pred.ben.rel.raw <- ben.obj$pred.ben.rel.raw
  risk.flipped.w <- ben.obj$risk.flipped.w
  risk.regular.w <- ben.obj$risk.regular.w
  
  # fit baseline model
  baseline.mod  <- baseline.risk(X = X, y = y, alpha = alpha)
  baseline.risk <- baseline.mod$response 
  
  # return
  return(list(
    inputs = list(X = X, w = w, y = y.orig, 
                  lifeyears = lifeyears, 
                  prediction.timeframe = prediction.timeframe, 
                  y.prediction.timeframe = y),
    models = list(glmnet.object = mod,
                  final.model = final.model,
                  coefficients.glmnet = structure(as.numeric(coefs.obj), names = c("(Intercept)", colnames(Z))),
                  coefficients.final = final.model$coefficients),
    average.treatment.effect = mean(pred.ben.abs),
    risk = list(risk.regular.w = risk.regular.w,
                risk.flipped.w = risk.flipped.w,
                risk.baseline = baseline.risk),
    benefits = list(predicted.absolute.benefit = pred.ben.abs,
                    predicted.relative.benefit = pred.ben.rel,
                    predicted.absolute.benefit.raw = pred.ben.abs.raw,
                    predicted.relative.benefit.raw = pred.ben.rel.raw),
    C.statistics = list(c.index.outcome = C.index.outcome(y = y, risk.prediction = baseline.risk),
                        c.index.benefit = C.index.benefit(y = y, w = w, predicted.benefit = pred.ben.abs))
  ))
  
} # FUN


effect.modeling_imputation.accounter <- function(predictive.model.imputed){
  
  # initialize
  pred.model.imp.adj <- list()
  m <- length(predictive.model.imputed)
  
  # ATE
  pred.model.imp.adj$average.treatment.effect <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$average.treatment.effect))
  
  # baseline model
  pred.model.imp.adj$baseline.model <- imputation.accounter_location(
    lapply(1:m, function(i){
      
      # get names of all variables (pre-selection)
      nam.all.variables <- c("(Intercept)", colnames(predictive.model.imputed[[i]]$inputs$X))
      
      # initialize long array with zeros for unselected variables
      selected.variables.long <- rep(0.0, length(nam.all.variables))
      names(selected.variables.long) <- nam.all.variables
      
      # assign values to long array
      selected.variables.short <- predictive.model.imputed[[i]]$baseline.model$coefficients
      selected.variables.long[names(selected.variables.short)] <- selected.variables.short
      selected.variables.long
      
    }))
  
  # effect model
  pred.model.imp.adj$effect.model <- imputation.accounter_location(
    lapply(1:m, function(i){
      
      # get names of all variables (pre-selection)
      nam.all.variables <- 
        names(predictive.model.imputed[[i]]$effect.model$model.selection$design.matrix_pre.selection)
      
      # initialize long array with zeros for unselected variables
      selected.variables.long <- rep(0.0, length(nam.all.variables))
      names(selected.variables.long) <- nam.all.variables
      
      # assign values to long array
      selected.variables.short <- predictive.model.imputed[[i]]$effect.model$summary[,1]
      selected.variables.long[names(selected.variables.short)] <- selected.variables.short
      selected.variables.long
      
    }))
  
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
  
  # C index outcome
  pred.model.imp.adj$C.statistics$c.index.outcome <- 
    imputation.accounter_scalar.location.stderr(lapply(1:m, function(i) predictive.model.imputed[[i]]$C.statistics$c.index.outcome))
  
  # C index benefit
  pred.model.imp.adj$C.statistics$c.index.benefit <- 
    imputation.accounter_scalar.location.stderr(lapply(1:m, function(i) predictive.model.imputed[[i]]$C.statistics$c.index.benefit))
  
  # return
  return(pred.model.imp.adj)
  
} # FUN
