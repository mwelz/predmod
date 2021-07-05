# load for baseline risk
source(paste0(getwd(), "/funs/linear-models/risk-modeling.R"))

# load imputation accounters
source(paste0(getwd(), "/funs/imputation/imputation.R"))


#' get the matrix "w * X[, interactions]"
#' 
#' @param X design matrix
#' @param w vector of binary treatment assignments
#' @param interacted.variables array of strings of the variables in X that shall be interacted with w. If equal to "all", all variables in X are interacted with _w_.
#' 
#' @export
get.interaction.terms.matrix <- function(X, w, interacted.variables){
  
  # error handling
  if(!is.character(interacted.variables)){
    stop("interacted.variables needs to be a vector of varable names of X.")
  }
  
  if("w" %in% interacted.variables) stop("w cannot be interacted with itself!")
  
  # add interaction variables
  interaction.terms.matrix <- sapply(which(colnames(X) %in% interacted.variables),
                                     function(j) ifelse(w == 1, X[,j], 0))
  
  colnames(interaction.terms.matrix) <- paste0("w.", interacted.variables)
  
  return(interaction.terms.matrix)
  
} # FUN


get.lambdapath <- function(y, x){
  
  # taken from https://stackoverflow.com/questions/23686067/default-lambda-sequence-in-glmnet-for-cross-validation
  
  n <- length(y)
  p <- ncol(x)
  
  ## Standardize variables: (need to use n instead of (n-1) as denominator)
  mysd <- function(z) sqrt(sum((z-mean(z))^2)/length(z))
  sx <- scale(x, scale = apply(x, 2, mysd))
  sx <- as.matrix(sx, ncol = n, nrow = p)
  
  ## Calculate lambda path (first get lambda_max):
  lambda_max <- max(abs(colSums(sx*y)))/n
  epsilon <- 0.0001
  K <- 100
  lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), 
                              length.out = K)), digits = 10)
  return(lambdapath)
} # FUN



#' helper function that recovers the index of the variables that shall always be retained
#' 
#' @param retained.variables A string array of variable names in Z that shall always be retained (also works on interaction variables). If NULL, then no restriction applies. Note that treatment assignment w will always be retained by the function.
#' @export
get.Z.index.of.retained.variables <- function(retained.variables, X, Z){
  
  Z.names <- colnames(Z)
  
  if(is.null(retained.variables)) return(1) # w will always be retained
  
  if(!any(retained.variables %in% Z.names)){
    stop("This variable/interaction effect does not exist in X")
  } 

  # w will always be retained, which is at index 1
  return(which(Z.names %in% c("w", retained.variables)))
  
} # FUN
  

#' Performs variable selection based on post-selection hypothesis tests. Function follows the strategy of Wasserman and Roeder (2009, Annals of Statistics):
#' 
#' step 0.1: create matrix Z = (w, X, interaction.terms). 
#' 
#' step 0.2: partition the observations in three roughly equally-sized sets sets: D1, D2, D3.
#' 
#' step 1.1: use D1 to perform variable selection by penalized logistic regression for given value of lambda. Let Z_lambda be the design matrix associated with the retained variables.
#' 
#' step 1.2: use the retained variables from 1.1 to calculate the logistic least squares estimator beta_lambda on D1.
#' 
#' step 1.3 calculate the empirical cross-entropy loss on D2, by using beta_lambda
#' 
#' step 2: repeat steps 1.1 to 1.3 for many lambdas. Let beta_lambdahat be the logistic least squares estimator that corresponds to the retained variables of the model with the lambda that minimizes the loss. Calculate this least squares estimator on D3.
#' 
#'  step 3: perform standard hypothesis tests to decide which of the retained variables in step 2 make it to the final model. Note that the critical value needs to be adjusted; see  Wasserman and Roeder (2009, Annals of Statistics) for details
#' 
#' @param X design matrix, can also be a data frame
#' @param y vector of binary responses. 
#' @param w vector of binary treatment assignments
#' @param interacted.variables string array of variables in _X_ that shall be interacted with treatment _w_
#' @param alpha the alpha as in 'glmnet'. Default is 1, which corresponds to the Lasso
#' @param retained.variables string array of variable names in Z that shall always be retained (i.e. also interaction terms can be considered). If NULL (default), then no restriction applies. Note that treatment assignment w will always be retained by the function.
#' @param significance.level for the hypothesis tests. Default is 0.05
#' 
#' @export
variable.selection <- function(X, y, w,
                               interacted.variables,
                               alpha = 1,
                               retained.variables = NULL,
                               significance.level = 0.05){
  
  ## prepare the design matrix Z 
  n    <- nrow(X)
  p    <- ncol(X)
  
  if(is.null(colnames(X))) colnames(X) <- paste0("X", 1:p)
  
  # get design matrix for the modeling
  if(is.null(interacted.variables)){
    Z <- data.frame(w = w, X)
  } else{
    Z <- data.frame(w = w, X,
                    get.interaction.terms.matrix(X = X, 
                                                 w = w, 
                                                 interacted.variables = interacted.variables))
  } # IF
  
  # must be a matrix
  Z <- as.matrix(Z)
  
  # get Z index of forcefully retained variables
  retained.variables.Z.idx <- get.Z.index.of.retained.variables(retained.variables = retained.variables, 
                                                                X = X, Z = Z)
  
  # randomly split data into three roughly equally sized samples
  set1 <- sample(1:n, floor(n / 3), replace = FALSE)
  set2 <- sample(setdiff(1:n, set1), floor(n / 3), replace = FALSE)
  set3 <- setdiff(1:n, c(set1, set2))
  
  # prepare cross-validation: loop over the lambdas of the glmnet path
  lambdapath <- get.lambdapath(x = Z[set1,], y = y[set1])
  
  suite <- lapply(1:length(lambdapath), function(...) list(retained.variables = NA, lasso.estimates = NA))
  loss  <- rep(NA_real_, length(lambdapath))
  
  for(i in 1:length(lambdapath)){
    
    # penalized logistic regression on first set
    mod <- glmnet::glmnet(Z[set1,], y[set1], 
                          family = "binomial", 
                          alpha = alpha,
                          lambda = lambdapath[i])
    
    # obtain retained variables (excluding intercept)
    kept.vars     <- glmnet::coef.glmnet(mod)@i 
    if(0 %in% kept.vars){
      kept.vars <- kept.vars[-which(kept.vars == 0)]
    } # IF
    
    # make sure that all variables that are forced to be retained will be retained
    kept.vars <- unique(sort(c(kept.vars, retained.variables.Z.idx), decreasing = FALSE))
    
    # add to suite
    suite[[i]]$retained.variables <- kept.vars
    suite[[i]]$lasso.estimates    <- c(as.numeric(mod$a0), as.numeric(mod$beta))
    
    #  calculate logistic least squares estimator with retained variables on first set
    df <- data.frame(y = y[set1], Z[set1, kept.vars])
    colnames(df) <- c("y", colnames(Z)[kept.vars])
    mod <- glm(y~., family =  binomial(link = "logit"), data = df)
    
    # predict Pr(Y=1) on second set
    df <- data.frame(Z[set2, kept.vars])
    colnames(df) <- colnames(Z)[kept.vars]
    p.hat <- predict.glm(mod, newdata = df, type = "response")
    
    # cross-validation: evaluate the loss (cross-entropy here) on the second set
    crossentropy <- rep(NA_real_, length(set2))
    crossentropy[y[set2] == 1] <- -log(p.hat[y[set2] == 1])
    crossentropy[y[set2] == 0] <- -log(1 - p.hat[y[set2] == 0])
    loss[i] <- mean(crossentropy)
    
  } # FOR
  
  # find the lambda that minimizes the empirical loss and its associated model
  lambda.min <- lambdapath[which.min(loss)]
  S.hat      <- suite[[which.min(loss)]]$retained.variables
  Z.lambda.min <- Z[, S.hat]
  
  # on third set: use S.hat to calculate logistic least squares estimator with Z.lambda.min
  df <- data.frame(y = y[set3], Z.lambda.min[set3,])
  colnames(df) <- c("y", colnames(Z.lambda.min))
  mod <- glm(y~., family =  binomial(link = "logit"), 
             data = df)
  
  # on third set: hypothesis testing
  coeffs <- summary(mod)$coefficients
  coeffs <- cbind(coeffs, rep(NA_real_, length(S.hat) + 1))
  colnames(coeffs) <- c("Estimate", "Std. Error", "z value", "critical value", "retain?")
  critval <- qnorm(significance.level / 2 * length(S.hat), lower.tail = FALSE)
  coeffs[, "critical value"] <- critval
  coeffs[, "retain?"] <- 1 * (abs(coeffs[, "z value"]) > critval)
  
  # final model
  D.hat   <- which(coeffs[-1, "retain?"] == 1)
  
  # make sure that all variables that are forced to be retained will be retained
  D.hat <- unique(sort(c(D.hat, retained.variables.Z.idx), decreasing = FALSE))
  
  # get the regularized estimates at the minimizing lambda
  regularized.estimates_lambda.min <- suite[[which.min(loss)]]$lasso.estimates
  names(regularized.estimates_lambda.min) <- c("(Intercept)", colnames(Z))
  
  # organize output
  model.selection <- list(selected.variables = colnames(Z)[D.hat],
                          design.matrix_selected.variables = Z[, D.hat],
                          formula.final.model = paste0("y ~ ", paste(colnames(Z)[D.hat], collapse = " + ")),
                          final.model.object = mod,
                          lambda.min = lambda.min,
                          significance.tests_lambda.min = coeffs,
                          regularized.estimates_lambda.min = regularized.estimates_lambda.min,
                          S.hat = S.hat)
  partitioning.membership <- list(set1 = sort(set1, decreasing = FALSE),
                                  set2 = sort(set2, decreasing = FALSE),
                                  set3 = sort(set3, decreasing = FALSE))
  
  return(list(final.model = model.selection, 
              partitioning.membership = partitioning.membership))
} # FUN



#' helper function that calculates the predicted benefits
#' 
#' @param X covariates
#' @param y binary outcomes
#' @param w binary treatment assignment
#' @param final.model final model after post-selection inference
#' @param Z covariates corresponding to the final model
#' 
#' @export
effect.model.predicted.benefits <- function(X, y, w, final.model, Z){
  
  # get risk with regular w (response is a probability)
  response <- unname(predict.glm(final.model,
                                 newdata = as.data.frame(Z), 
                                 type = "response"))
  
  # get risk with flipped w
  kept.interactions <- substr(colnames(Z)[startsWith(colnames(Z), "w.")], start = 3, stop = 1e6)
  kept.X            <- colnames(Z)[!startsWith(colnames(Z), "w.")]
  kept.X            <- kept.X[-which(kept.X == "w")] # w doesn't belong to X
  interaction.terms.flipped <- sapply(kept.interactions, 
                                      function(j) ifelse(w == 1, 0, X[,j]))
  w.flipped   <- ifelse(w == 1, 0, 1)
  
  if(length(kept.interactions) == 0){
    Z_w.flipped <- data.frame(w = w.flipped, X[,kept.X])
  } else{
    Z_w.flipped <- data.frame(w = w.flipped, X[,kept.X], interaction.terms.flipped)
  } # IF
  colnames(Z_w.flipped) <- colnames(Z)
  
  response_w.flipped <- unname(predict.glm(final.model, 
                                           newdata = as.data.frame(Z_w.flipped),
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


#' Performs effect modeling. Final effect model is done via variable selection based on post-selection hypothesis tests. Function follows the strategy of Wasserman and Roeder (2009, Annals of Statistics):
#' 
#' step 0.1: create matrix Z = (w, X, interaction.terms). 
#' 
#' step 0.2: partition the observations in three roughly equally-sized sets sets: D1, D2, D3.
#' 
#' step 1.1: use D1 to perform variable selection by penalized logistic regression for given value of lambda. Let Z_lambda be the design matrix associated with the retained variables.
#' 
#' step 1.2: use the retained variables from 1.1 to calculate the logistic least squares estimator beta_lambda on D1.
#' 
#' step 1.3 calculate the empirical cross-entropy loss on D2, by using beta_lambda
#' 
#' step 2: repeat steps 1.1 to 1.3 for many lambdas. Let beta_lambdahat be the logistic least squares estimator that corresponds to the retained variables of the model with the lambda that minimizes the loss. Calculate this least squares estimator on D3.
#' 
#'  step 3: perform standard hypothesis tests to decide which of the retained variables in step 2 make it to the final model. Note that the critical value needs to be adjusted; see  Wasserman and Roeder (2009, Annals of Statistics) for details
#' 
#' @param X design matrix, can also be a data frame
#' @param y vector of binary responses. 
#' @param w vector of binary treatment assignments
#' @param interacted.variables string array of variables in _X_ that shall be interacted with treatment _w_
#' @param alpha the alpha as in 'glmnet'. Default is 1, which corresponds to the Lasso
#' @param lifeyears vector of life years. Default is NULL.
#' @param prediction.timeframe vector of the prediction time frame. Default is NULL.
#' @param retained.variables string array of variable names in Z that shall always be retained (i.e. also interaction terms can be considered). If NULL (default), then no restriction applies. Note that treatment assignment w will always be retained by the function.
#' @param significance.level for the hypothesis tests. Default is 0.05
#' 
#' @export
effect.modeling <- function(X, y, w,
                            interacted.variables,
                            alpha = 1,
                            lifeyears = NULL, 
                            prediction.timeframe = NULL,
                            retained.variables = NULL,
                            significance.level = 0.05){
  
  # truncate y if necessary
  y.orig    <- y
  
  if(!is.null(lifeyears) & !is.null(prediction.timeframe)){
    lifeyears <- ifelse(lifeyears <= prediction.timeframe, lifeyears, prediction.timeframe) 
    y         <- ifelse(lifeyears <= prediction.timeframe, y, 0)
  } # IF
  
  ### 1. fit baseline risk  ----
  # no information on w allowed, so we cannot use the retained variables from the effect modeling
  baseline.mod  <- baseline.risk(X = X, y = y, alpha = 1)
  baseline.risk <- baseline.mod$response 
  
  
  ### 2. perform variable selection ----
  # We use the strategy in Wasserman and Roeder (2009; Annals of Statistics)
  vs.obj <- variable.selection(X = X, y = y, w = w, 
                               interacted.variables = interacted.variables, 
                               alpha = alpha, 
                               retained.variables = retained.variables, 
                               significance.level = significance.level)
  
  # get design matrix associated with the final model
  Z     <- vs.obj$final.model$design.matrix_selected.variables
  
  ### 3. fit the final effect model and use it for risk estimates ----
  # fit the final model, this time on full sample (no penalty required, selection already took place!)
  final.model        <- glm(y ~., data = data.frame(y, Z), 
                            family =  binomial(link = "logit"))
  predicted.benefits <- effect.model.predicted.benefits(X = X, y = y, w = w, 
                                                        final.model = final.model, Z = Z)
  
  # coefficients
  coefficients <- summary(final.model)$coefficients
  
  
  ### 4. calculate the C statistics ----
  
  # match cases based on observed benefit TODO: fix the data thing everywhere!
  matched <- MatchIt::matchit(w ~ pb, data = data.frame(w=w, pb=predicted.benefits$pred.ben.abs))
  match.treated <- as.numeric(rownames(matched$match.matrix))
  match.control <- as.numeric(matched$match.matrix[,1])
  
  # remove unpaired observations
  no.pairing <- which(is.na(match.control))
  if(length(no.pairing) > 0){
    match.treated <- match.treated[-no.pairing]
    match.control <- match.control[-no.pairing]
  }
  
  
  # observed benefit & C index for benefit
  obs.ben             <- y[match.control] - y[match.treated]
  pred.ben.abs.paired <- (predicted.benefits$pred.ben.abs[match.control] +
                            predicted.benefits$pred.ben.abs[match.treated]) / 2
  c.index.benefit     <- unname(Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)[1])
  
  # C index
  c.index.outcome <- unname(Hmisc::rcorr.cens(predicted.benefits$risk.regular.w, y)[1])
  
  
  ### 5. return ----
  return(list(
    inputs = list(X = X, w = w, y = y.orig, 
                  lifeyears = lifeyears, 
                  prediction.timeframe = prediction.timeframe, 
                  y.prediction.timeframe = y),
    average.treatment.effect = mean(predicted.benefits$pred.ben.abs),
    baseline.model = baseline.mod,
    effect.model = list(formula = vs.obj$final.model$formula.final.model, 
                        selected.data = Z,
                        glm.obj = final.model,
                        coefficients = coefficients[,1],
                        summary = coefficients,
                        model.building = vs.obj),
    risk = list(risk.regular.w = predicted.benefits$risk.regular.w,
                risk.flipped.w = predicted.benefits$risk.flipped.w,
                risk.baseline  = baseline.risk),
    benefits = list(predicted.absolute.benefit = predicted.benefits$pred.ben.abs,
                    predicted.relative.benefit = predicted.benefits$pred.ben.rel,
                    predicted.absolute.benefit.raw = predicted.benefits$pred.ben.abs.raw,
                    predicted.relative.benefit.raw = predicted.benefits$pred.ben.rel.raw),
    C.statistics = list(c.index.outcome = c.index.outcome,
                        c.index.benefit = c.index.benefit)
  ))
} # FUN


#' TODO: write documentation and come up with way to meaningfully return the effect models
#'
#'
effect.modeling_imputation.accounter <- function(predictive.model.imputed){
  
  # initialize
  pred.model.imp.adj <- list()
  m <- length(predictive.model.imputed)
  
  # ATE
  pred.model.imp.adj$average.treatment.effect <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$average.treatment.effect))
  
  # baseline model
  pred.model.imp.adj$baseline.model <- "TODO"
  
  # effect model
  pred.model.imp.adj$effect.model <- "TODO"
  
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
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$C.statistics$c.index.outcome))
  
  # C index benefit
  pred.model.imp.adj$C.statistics$c.index.benefit <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$C.statistics$c.index.benefit))
  
  # return
  return(pred.model.imp.adj)
  
} # FUN