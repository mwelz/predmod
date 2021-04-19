cox.effect.modeling <- function(X, w, y, lifeyears, predictiontimeframe,
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
  
  # add interaction variables
  interaction.terms <- sapply(which(colnames(X) %in% interaction.vars),
                              function(j) ifelse(w == 1, X[,j], 0))
  
  interaction.terms <- cbind(w, interaction.terms)
  colnames(interaction.terms) <- c("w", paste0("w.", interaction.vars))
  X.star <- cbind(X, interaction.terms)
  glmsurvival.coxeffect.obj <-survival::Surv(lifeyears, y)
  colnames(glmsurvival.coxeffect.obj) <- c("time", "status")
  #cox.models.stage1 <- glmnet::cv.glmnet(x, glmsurvival.coxeffect.obj, family = "cox", type.measure = "C", alpha = 1)
  
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
  final.model.fit <- rms::cph(glmsurvival.coxeffect.obj~ X.final,y=TRUE,x=TRUE)
  lp=final.model.fit$linear.predictors
  basehazard <- survival::basehaz(final.model.fit,centered=FALSE)
  Hazard_t_index = match(1, round(basehazard$time,1) == predictiontimeframe, nomatch = NA)
  lp2<- as.numeric(X.final  %*% final.model.fit$coefficients)
  #TO CHECK: lp2 needs to be equal to lp; but there seems to be a mismatch
  # 1) ordinary w
  probs  =  1-exp(-basehazard[Hazard_t_index,1])^exp(lp) 

  
  
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
}
