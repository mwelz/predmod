

effect.modeling.poiss <- function(X, w, y,ly, 
                                  alpha = alpha, 
                                  interactions = NULL,
                                  sig.level = 0.05, ...){
  
  ### 0. preparation ----
  # split the sample as suggested in Wasserman and Roeder (2009)
  n    <- nrow(X)
  p    <- ncol(X)
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
  
  #Set X to X-star
  X.star <- X
  
  # add interaction variables
  if(is.null(interactions)==TRUE || is.na(interactions)==FALSE)
  {
    interaction.terms <- sapply(which(colnames(X) %in% interaction.vars),
                                function(j) ifelse(w == 1, X[,j], 0))
    
    interaction.terms <- cbind(w, interaction.terms)
    colnames(interaction.terms) <- c("w", paste0("w.", interaction.vars))
    X.star <- cbind(X, interaction.terms)
  }
  ### 1. stage 1: penalized regression on whole set  ----
  # whole set is x.star. We apply sample splitting as suggested in Wasserman and Roeder (2009)
  mod.pm <- glmnet::cv.glmnet(X.star[set1,], y[set1], family = "poisson",offset=log(ly[set1]), alpha = alpha)
  
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
  fit.1     <- glm(y ~., data = data.frame(y, X.kept)[set2,], family = "poisson",offset=log(ly[set2]))
  fit.0     <- glm(y ~., data = data.frame(y, X.star.red)[set2,], family = "poisson",offset=log(ly[set2])) # reduced model
  formula.0 <- paste0("y ~ ", paste(colnames(X.star.red), collapse = " + "))
  formula.1 <- paste0("y ~ ", paste(colnames(X.kept), collapse = " + "))
  
  # perform likelihood ratio test
  test.stat <- fit.0$deviance - fit.1$deviance
  df        <- fit.0$df.residual - fit.1$df.residual
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
  probs <- unname(predict.glm(final.model,
                              newdata = as.data.frame(X.final), 
                              type = "response"))
  
  # get risk with flipped w
  kept.X <- colnames(X.final)[startsWith(colnames(X.final), "X")]
  kept.w <- sub(".*\\.", "", colnames(X.final)[startsWith(colnames(X.final), "w")])
  X.star <- cbind(X, w = ifelse(w == 1, 0, 1))
  interaction.terms.flipped <- sapply(kept.w, 
                                      function(j) ifelse(w == 1, 0, X.star[,j]))
  X.rev <- cbind(X[,kept.X], interaction.terms.flipped)
  colnames(X.rev) <- colnames(X.final)
  probs.flipped.w <- unname(predict.glm(final.model, 
                                        newdata = as.data.frame(X.rev),
                                        type = "response"))
  
  # get observed predicted benefit
  pred.ben.abs.raw <- probs - probs.flipped.w
  pred.ben.abs     <- ifelse(w == 1, -pred.ben.abs.raw, pred.ben.abs.raw)
  Estimated_reduction_effect_per_1000 =  (sum(pred.ben.abs_effect)/sum(ly))*1000
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
  formula.final.model                   <- paste0("y ~ ", 
                                                  paste(colnames.final, collapse = " + "))
  
  # coefficients
  coefficients <- final.model$coefficients
  names(coefficients) <- c("(Intercept)", colnames.final)
  
  # summary
  summry <- summary(final.model)$coefficients
  rownames(summry) <- names(coefficients)
  
  ## 5. fit baseline risk  ----
  # no information on w allowed, so we cannot use the retained variables from the effect modeling
  baseline.mod <- risk.model.poiss.stage1(X = X, y = y, ly=ly, alpha = alpha)
  basline.poiss <-   baseline.mod$lp
  
  # calculate C index by using predicted risk (with regular w)
  c.index <- unname(Hmisc::rcorr.cens(x = probs, S = y)[1])
  
  ## 6. return ----
  return(list(
    inputs = list(X = X, w = w, y = y,ly=ly),
    baseline.model = baseline.mod,
    effect.model = list(formula = formula.final.model, 
                        selected.data = X.final,
                        glm.obj = final.model,
                        coefficients = coefficients,
                        summary = summry,
                        model.building = list(test.stat = test.stat,
                                              df = df,
                                              pval = pval)),
    risk.regular.w = probs,
    risk.flipped.w = probs.flipped.w,
    basline.poiss = basline.poiss,
    predicted.absolute.benefit = pred.ben.abs,
    predicted.absolute.benefit.raw = pred.ben.abs.raw,
    predicted.relative.benefit = pred.ben.rel,
    predicted.relative.benefit.raw = pred.ben.rel.raw,
    ate.hat = mean(pred.ben.abs),
    ratereduction.1000ly = Estimated_reduction_per_1000_lY,
    c.index =NA
    # c.index = c.index
  ))
}

