source(paste0(getwd(), "/funs/c-statistics/c-statistics.R"))
source(paste0(getwd(), "/funs/linear-models/effect-modeling.R"))
source(paste0(getwd(), "/funs/cox/cox-risk-modeling.R"))


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
#' @param prediction.timeframe
#' @param lifeyears
#' 
#' @export
variable.selection.cox <- function(X, y, w,
                                   interacted.variables,
                                   alpha = 1,
                                   prediction.timeframe,
                                   lifeyears,
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
    
    # get survival matrix; to be passed as response of glmnet
    survival.matrix           <-survival::Surv(lifeyears[set1], y[set1])
    colnames(survival.matrix) <- c("time", "status")
    
    # penalized logistic regression on first set
    mod <- glmnet::glmnet(x = Z[set1,], y = survival.matrix, 
                          family = "cox", type.measure = "C", 
                          alpha = alpha,
                          lambda = lambdapath[i])
    
    # obtain retained variables. Cox doesn't have intercept, so account for zero indexing
    kept.vars     <- glmnet::coef.glmnet(mod)@i + 1
    
    # make sure that all variables that are forced to be retained will be retained
    kept.vars <- unique(sort(c(kept.vars, retained.variables.Z.idx), decreasing = FALSE))
    
    # add to suite
    suite[[i]]$retained.variables <- kept.vars
    suite[[i]]$lasso.estimates    <- as.numeric(mod$beta) # Cox doesn't have intercept
    
    #  fit Cox model with retained variables on first set
    mod <- survival::coxph(survival.matrix~., data = data.frame(Z[set1, kept.vars]))
    # CONTINUE HERE
    
    
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



cox.effect.modeling <- function(X, y, w,
                                interacted.variables,
                                alpha = 1,
                                lifeyears, 
                                prediction.timeframe,
                                retained.variables = NULL,
                                significance.level = 0.05){
  
  # truncate y if necessary
  y.orig    <- y
  lifeyears <- ifelse(lifeyears <= prediction.timeframe, lifeyears, prediction.timeframe) 
  y         <- ifelse(lifeyears <= prediction.timeframe, y, 0)
  
  # make X a matrix
  X <- as.matrix(X)
  
  ### 1. fit baseline risk via Cox modeling  ----
  # no information on w allowed, so we cannot use the retained variables from the effect modeling
  baseline.mod  <- baseline.risk.cox(X = X, y = y, alpha = alpha, 
                                     lifeyears = lifeyears,
                                     prediction.timeframe = prediction.timeframe)
  baseline.risk <- baseline.mod$response 


  
  ###################################################################
  ### 0. preparation ----
  # split the sample as suggested in Wasserman and Roeder (2009)
  lifeyears <- ifelse(lifeyears <=predictiontimeframe, lifeyears, predictiontimeframe)
  y <- ifelse(lifeyears <=predictiontimeframe, y, 0)
  n    <- nrow(X)
  p    <- ncol(X)
  set.seed(25)
  set1 <- sample(1:n, floor(0.5 * n), replace = FALSE)
  set2 <- setdiff(1:n, set1)
  
  # prepare the variable set as in rekkas2019:
  if(!is.null(colnames(X))){
    colnames.orig <- colnames(X)
  } else{
    colnames.orig <- paste0("X", 1:p)
  }
  colnames(X)     <- paste0("X", 1:p)
  
  if(is.null(interactions)){
    interaction.vars <- colnames(X)
  } else{
    interaction.vars <- paste0("X", which(colnames.orig %in% interactions))
  }
  
  # add interaction variables
  interaction.terms <- sapply(which(colnames(X) %in% interaction.vars),
                              function(j) ifelse(w == 1, X[,j], 0))
  
  interaction.terms <- cbind(w, interaction.terms)
  colnames(interaction.terms) <- c("w", paste0("w.", interaction.vars))
  X.star <- cbind(X, interaction.terms)
  glmsurvival.coxeffect.obj <-survival::Surv(lifeyears, y)
  colnames(glmsurvival.coxeffect.obj) <- c("time", "status")
  
  ### 1. stage 1: penalized regression on whole set  ----
  # whole set is x.star. We apply sample splitting as suggested in Wasserman and Roeder (2009)
  mod.pm <- glmnet::cv.glmnet(X.star[set1,], glmsurvival.coxeffect.obj[set1], family = "cox", type.measure = "C", alpha = 1)
  
  ### 2. stage 2: perform variable selection based on the "best" model ----
  coefs.obj     <- glmnet::coef.glmnet(mod.pm, s = "lambda.min")
  kept.vars     <- coefs.obj@i 
  if(0 %in% kept.vars){
    kept.vars <- kept.vars[-which(kept.vars == 0)]
  }
  kept.vars.nam <- colnames(X.star)[kept.vars]
  
  ### 3. stage 3: perform chi-squared test on significance of interactions in the "best" model ----
  # big model
  X.kept <- X.star[, kept.vars.nam]
  
  # reduced model: (X_lambda, w), where X_lambda are the retained (by the lasso) initial variables
  kept.vars.no.x.nam <- 
    kept.vars.nam[grepl(pattern = "w.", x = kept.vars.nam, fixed = TRUE)]
  X.star.red <- X.kept[, -which(kept.vars.nam %in% kept.vars.no.x.nam)]
  
  # fit both large and reduced model
  #fit.1     <-  glmnet::glmnet(X.kept[set2,], glmsurvival.coxeffect.obj[set2], family = "cox", type.measure = "C")
  #fit.0     <-  glmnet::glmnet(X.star.red[set2,], glmsurvival.coxeffect.obj[set2], family = "cox", type.measure = "C")
  fit.1     <- rms::cph(glmsurvival.coxeffect.obj[set2]~ X.kept[set2,])
  fit.0     <-  rms::cph(glmsurvival.coxeffect.obj[set2]~ X.star.red[set2,])# reduced model
  formula.0 <- paste0("surv ~ ", paste(colnames(X.star.red), collapse = " + "))
  formula.1 <- paste0("surv ~ ", paste(colnames(X.kept), collapse = " + "))
  anova_fit1 = anova(fit.1)
  anova_fit0 = anova(fit.0)
  
  # perform likelihood ratio test
  test.stat <- anova_fit1[2,1]- anova_fit0[2,1]
  df        <- anova_fit1[2,2]- anova_fit0[2,2]
  pval      <- pchisq(test.stat, df, lower.tail = FALSE)
  
  # if we reject the null, go with the smaller (the reduced model) to make risk predictions
  if(pval < sig.level){
    final.model         <- fit.0
    X.final             <- X.star.red
  } else{
    final.model         <- fit.1
    X.final             <- X.kept
  }
  
  ### 4. use the final model to obtain risk estimates ----
  
  
  # risk with regular w
  #final.model.fit <-glmnet::glmnet(X.final, glmsurvival.coxeffect.obj, family = "cox", type.measure = "C")
  #lambda        <- min(final.model.fit$lambda)
  #lp              <-predict(final.model.fit, s=lambda, X.final,type="link") #linear predictors
  #basehaz  <- hdnom::glmnet_basesurv(lifeyears, y, lp , times.eval = predictiontimeframe, centered = FALSE)
  final.model.fit <- rms::cph(glmsurvival.coxeffect.obj~ X.final,y=TRUE,x=TRUE)
  lp=X.final %*% final.model.fit$coefficients
  basehazard <- survival::basehaz(final.model.fit,centered=FALSE)
  Hazard_t_index = match(1, round(basehazard$time,1) == predictiontimeframe, nomatch = NA)
  
  # 1) ordinary w
  probs  =  1-exp(-basehazard[Hazard_t_index,1])^exp(lp) 
  #probs  =  1-exp(-basehaz$cumulative_base_hazard)^exp(lp)
  # coefs.glm.cox     <- glmnet::coef.glmnet(final.model.fit, s = lambda)
  
  # get risk with flipped w
  kept.X <- colnames(X.final)[startsWith(colnames(X.final), "X")]
  kept.w <- sub(".*\\.", "", colnames(X.final)[startsWith(colnames(X.final), "w")])
  X.star <- cbind(X, w = ifelse(w == 1, 0, 1))
  interaction.terms.flipped <- sapply(kept.w, 
                                      function(j) ifelse(w == 1, 0, X.star[,j]))
  X.rev <- cbind(X[,kept.X], interaction.terms.flipped)
  colnames(X.rev) <- colnames(X.final)
  
  lp.flipped.w     <- as.numeric(X.rev  %*% final.model.fit$coefficients)
  probs.flipped.w  <-  1-exp(-basehazard[Hazard_t_index,1])^exp(lp.flipped.w) 
  
  
  # get absolute predicted benefit
  pred.ben.abs.raw <- probs - probs.flipped.w
  pred.ben.abs     <- ifelse(w == 1, pred.ben.abs.raw, -pred.ben.abs.raw)
  
  # get relative predicted benefit
  pred.ben.rel.raw <- probs / probs.flipped.w
  pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)
  
  
  ### 5. housekeeping: make sure the naming is consistent
  # update variable names
  colnames.final.temp <- colnames(X.final)
  colnames.final      <- rep(NA_character_, length(colnames.final.temp))
  
  for(i in 1:p){
    idx <- grepl(pattern = paste0("X", i), x = colnames.final.temp, fixed = TRUE)
    colnames.final[idx] <- gsub(paste0("X", i), colnames.orig[i], colnames.final.temp[idx])
  }
  colnames.final[is.na(colnames.final)] <- "w" # fill up 'w'
  colnames(X.final)                     <- colnames.final
  colnames(X)                           <- colnames.orig
  formula.final.model                   <- paste0("surv ~ ", 
                                                  paste(colnames.final, collapse = " + "))
  
  # coefficients
  coefficients <- final.model$coefficients
  names(coefficients) <- c(colnames.final)
  
  
  
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
  Dxy = final.model$stats[[9]]
  abs.Dxy = abs(Dxy)
  c.index.youtcome=(abs.Dxy/2)+0.5
  
  
  ## 6. return ----
  return(list(
    inputs = list(X = X, w = w, y = y),
    effect.model = list(formula = formula.final.model, 
                        selected.data = X.final,
                        glm.obj = final.model,
                        coefficients = coefficients,
                        model.building = list(test.stat = test.stat,
                                              df = df,
                                              pval = pval)),
    risk.regular.w = probs,
    risk.flipped.w = probs.flipped.w,
    predicted.absolute.benefit = pred.ben.abs,
    predicted.absolute.benefit.raw = pred.ben.abs.raw,
    predicted.relative.benefit = pred.ben.rel,
    predicted.relative.benefit.raw = pred.ben.rel.raw,
    ate.hat = mean(pred.ben.abs),
    c.index.youtcome = c.index.youtcome,
    c.index.benefit =c.index.benefit
  ))
  
  
} # FUN