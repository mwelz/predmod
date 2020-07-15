
transform.to.probability <- function(x) 1 / (1 + exp(-x))


quantile.group <- function(x, cutoffs = c(0.25, 0.5, 0.75)){
  # cutoffs are the quantile cutoffs (like c(0.25, 0.5, 0.75))
  q         <- quantile(x, cutoffs)
  q         <- c(-Inf, q, Inf)
  groups    <- as.character(cut(x, breaks = q, include.lowest = FALSE, right = TRUE, dig.lab = 3))
  group.nam <- sort(unique(groups), decreasing = TRUE)
  group.mat <- matrix(NA, length(x), length(group.nam))
  nam       <- rep(NA, length(group.nam))
  
  for(j in 1:length(group.nam)){
    if(j == 1){
      nam[j] <- paste0("<=", 100*cutoffs[j], "% quantile")
    } else if (j == length(group.nam)){
      nam[j] <- paste0(">", 100*cutoffs[j-1], "% quantile")
    } else{
      nam[j] <- paste0("(", 100*cutoffs[j-1], ",", 100*cutoffs[j], "]% quantile")
    }
    group.mat[,j] <- groups == group.nam[j]
  }
  colnames(group.mat) <- nam
  return(group.mat)
}

risk.model.stage1 <- function(X, y, alpha = 1){
  
  # stage 1 doesn't use W!
  models.stage1 <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = alpha)
  lambda        <- models.stage1$lambda.min # best lambda (needs to be fixed in second stage)
  coefs.obj     <- glmnet::coef.glmnet(models.stage1, s = "lambda.min")
  kept.vars     <- coefs.obj@i
  if(0 %in% kept.vars) kept.vars <- kept.vars[-which(kept.vars == 0)]
  coefs         <- coefs.obj@x
  X.star        <- cbind(intercept = 1, X[,kept.vars])
  lp            <- as.numeric(X.star %*% coefs)
  
  return(list(
    models.stage1.cv = models.stage1,
    lambda.min = lambda,
    lp = lp
  ))
}


risk.model.stage2 <- function(lp, y, w, lambda, offset.lp = TRUE ){

  X.stage2 <- cbind(w = w, wlp = w * lp)
  
  if(offset.lp){
    oset <- lp
  } else{
    oset <- rep(0, length(lp))
  }
  
  # use same lambda as model 1!
  mod.stage2 <- glmnet::glmnet(X.stage2, y, 
                               family = "binomial",
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
  # 2) w flipped
  w.rev           <- ifelse(w == 1, 0, 1)
  X.stage2.rev    <- cbind(w = w.rev, wlp = w.rev * lp)
  probs.flipped.w <- glmnet::predict.glmnet(mod.stage2, 
                                            newx = X.stage2.rev,
                                            type = "response", 
                                            newoffset = oset)
  
  # return
  return(list(mod.stage2 = mod.stage2,
              risk.pred.regular = transform.to.probability(as.numeric(probs.regular)),
              risk.pred.flipped.w = transform.to.probability(as.numeric(probs.flipped.w))
  ))
}


risk.modeling <- function(X, w, y, alpha, offset.lp = TRUE){
  ## stage 1
  stage1 <- risk.model.stage1(X = X, y = y, alpha = alpha)
  
  ## stage 2
  stage2 <- risk.model.stage2(lp = stage1$lp,
                              y = y, 
                              w = w,
                              lambda = stage1$lambda.min,
                              offset.lp = offset.lp)
  
  # predicted benefit
  pred.ben.raw <- stage2$risk.pred.regular - stage2$risk.pred.flipped.w
  pred.ben     <- ifelse(w == 1, -pred.ben.raw, pred.ben.raw)
  
  # coefficients
  coefs.stage1 <- as.matrix(glmnet::coef.glmnet(stage1$models.stage1.cv, s = "lambda.min"))
  coefs.stage2 <- as.matrix(glmnet::coef.glmnet(stage2$mod.stage2))
  colnames(coefs.stage2) <- colnames(coefs.stage1) <- "Estimated Coefficient"
  
  
  return(list(
    inputs = list(X = X, w = w, y = y),
    mod.stage1 = stage1$models.stage1.cv,
    mod.stage2 = stage2$mod.stage2,
    coefficients.stage1 = coefs.stage1,
    coefficients.stage2 = coefs.stage2,
    linear.predictor = stage1$lp,
    risk.baseline = transform.to.probability(stage1$lp),
    risk.regular.w = stage2$risk.pred.regular,
    risk.flipped.w = stage2$risk.pred.flipped.w,
    predicted.benefit = pred.ben,
    predicted.benefit.raw = pred.ben.raw,
    ate.hat = mean(pred.ben)
  ))
}


get.benefits <- function(risk.model.obj, cutoffs = c(0.25, 0.5, 0.75)){
  y <- risk.model.obj$inputs$y
  w <- risk.model.obj$inputs$w
  
  # group observations by their quantile of predicted baseline risk
  quantile.groups <- quantile.group(risk.model.obj$risk.baseline, cutoffs)
  
  # extract predicted benefit
  pred.ben <- risk.model.obj$predicted.benefit
  
  ## calculate observed benefit and predicted benefit for each quantile group
  # initialize
  obs.ben.mat           <- as.data.frame(matrix(NA_real_, ncol(quantile.groups), 4))
  colnames(obs.ben.mat) <- c("quantile", "mean", "stderr", "df")
  obs.ben.mat$quantile  <- gsub(" .*$", "", colnames(quantile.groups))
  pred.ben.mat          <- obs.ben.mat
  
  for(i in 1:ncol(quantile.groups)){
    group <- quantile.groups[,i]
    
    ## observed benefit
    # corresponds to difference in mean(y[group & W=w])
    x1                        <- y[group & w == 1]
    x2                        <- y[group & w == 0]
    ttest                     <- t.test(x1, x2)
    obs.ben.mat[i, "mean"]    <- unname(ttest$estimate[1] - ttest$estimate[2])
    obs.ben.mat[i, "stderr"]  <- unname(ttest$stderr)
    obs.ben.mat[i, "df"]      <- unname(ttest$parameter)
    
    ## group by predicted benefit
    pred.ben.mat[i, "mean"]   <- mean(pred.ben[group])
    pred.ben.mat[i, "stderr"] <- sqrt(var(pred.ben[group]) / sum(group))
    pred.ben.mat[i, "df"]     <- sum(group) - 1
  }
  
  return(list(
    group.observed.benefit = obs.ben.mat,
    group.predicted.benefit = pred.ben.mat,
    group.membership = quantile.groups
  ))
}


calibration.plot <- function(risk.model.obj, quantiles = c(0.25, 0.5, 0.75), alpha.significance = 0.05){
  
  # get observed and predicted benefit by quantile group
  benefits <- get.benefits(risk.model.obj, cutoffs = quantiles)
  
  # make everything positive for visualization
  benefits$group.predicted.benefit$mean <- abs(benefits$group.predicted.benefit$mean) 
  benefits$group.observed.benefit$mean  <- abs(benefits$group.observed.benefit$mean)
  
  # limits of the plot
  limits <- c(min(benefits$group.predicted.benefit$mean, benefits$group.observed.benefit$mean, -0.2),
              max(benefits$group.predicted.benefit$mean, benefits$group.observed.benefit$mean, 0.3))
  
  # get the whiskers (SE*quantile)
  whisker.obs.ben <- benefits$group.observed.benefit$stderr *
    qt(1-alpha.significance/2, benefits$group.observed.benefit$df)
  
  # the plot
  library(ggplot2)
  ggplot(mapping = aes(x = benefits$group.predicted.benefit$mean,
                       y = benefits$group.observed.benefit$mean)) +
    geom_point() +
    geom_linerange(mapping = aes(ymin = benefits$group.observed.benefit$mean - whisker.obs.ben,
                                 ymax = benefits$group.observed.benefit$mean + whisker.obs.ben)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = limits, ylim = limits) +
    labs(x = "Predicted benefit", y = "Observed benefit") +
    theme_classic()
}


effect.modeling <- function(x, w, y, 
                            alpha = alpha, 
                            interactions = NULL,
                            sig.level = 0.05, ...){
  
  ### 0. preparation ----
  # split the sample as suggested in Wasserman and Roeder (2009)
  n    <- nrow(x)
  p    <- ncol(x)
  set1 <- sample(1:n, floor(0.5 * n), replace = FALSE)
  set2 <- setdiff(1:n, set1)
  
  # prepare the variable set as in rekkas2019:
  if(!is.null(colnames(x))){
    colnames.orig <- colnames(x)
  } else{
    colnames.orig <- paste0("x", 1:p)
  }
  colnames(x)     <- paste0("x", 1:p)
  
  if(is.null(interactions)){
    interaction.vars <- colnames(x)
  } else{
    interaction.vars <- paste0("x", which(colnames.orig %in% interactions))
  }
  
  # add interaction variables
  interaction.terms <- sapply(which(colnames(x) %in% interaction.vars),
                              function(j) ifelse(w == 1, x[,j], 0))
  
  interaction.terms <- cbind(w, interaction.terms)
  colnames(interaction.terms) <- c("w", paste0("w.", interaction.vars))
  x.star <- cbind(x, interaction.terms)
  
  ### 1. stage 1: penalized regression on whole set  ----
  # whole set is x.star. We apply sample splitting as suggested in Wasserman and Roeder (2009)
  mod.pm <- glmnet::cv.glmnet(x.star[set1,], y[set1], family = "binomial", alpha = alpha, ...)
  
  ### 2. stage 2: perform variable selection based on the "best" model ----
  coefs.obj     <- glmnet::coef.glmnet(mod.pm, s = "lambda.min")
  kept.vars     <- coefs.obj@i 
  if(0 %in% kept.vars){
    kept.vars <- kept.vars[-which(kept.vars == 0)]
  }
  kept.vars.nam <- colnames(x.star)[kept.vars]
  
  ### 3. stage 3: perform chi-squared test on significance of interactions in the "best" model ----
  # big model
  x.kept <- x.star[, kept.vars.nam]
  
  # reduced model: (X_lambda, w), where X_lambda are the retained (by the lasso) initial variables
  kept.vars.no.x.nam <- 
    kept.vars.nam[grepl(pattern = "w.", x = kept.vars.nam, fixed = TRUE)]
  x.star.red <- x.kept[, -which(kept.vars.nam %in% kept.vars.no.x.nam)]
  
  # fit both large and reduced model
  fit.1     <- glm(y ~., data = data.frame(y, x.kept)[set2,], family = binomial)
  fit.0     <- glm(y ~., data = data.frame(y, x.star.red)[set2,], family = binomial) # reduced model
  formula.0 <- paste0("y ~ ", paste(colnames(x.star.red), collapse = " + "))
  formula.1 <- paste0("y ~ ", paste(colnames(x.kept), collapse = " + "))
  
  # perform likelihood ratio test
  test.stat <- fit.0$deviance - fit.1$deviance
  df        <- fit.0$df.residual - fit.1$df.residual
  pval      <- pchisq(test.stat, df, lower.tail = FALSE)
  
  # if we reject the null, go with the smaller (the reduced model) to make risk predictions
  if(pval < sig.level){
    final.model         <- fit.0
    x.final             <- x.star.red
  } else{
    final.model         <- fit.1
    x.final             <- x.kept
  }
  
  ### 4. use the final model to obtain risk estimates ----
  
  # risk with regular w
  probs <- unname(predict.glm(final.model,
                              newdata = as.data.frame(x.final), 
                              type = "response"))
  
  # get risk with flipped w
  kept.x <- colnames(x.final)[startsWith(colnames(x.final), "x")]
  kept.w <- sub(".*\\.", "", colnames(x.final)[startsWith(colnames(x.final), "w")])
  x.star <- cbind(x, w = ifelse(w == 1, 0, 1))
  interaction.terms.flipped <- sapply(kept.w, 
                                      function(j) ifelse(w == 1, 0, x.star[,j]))
  x.rev <- cbind(x[,kept.x], interaction.terms.flipped)
  colnames(x.rev) <- colnames(x.final)
  probs.flipped.w <- unname(predict.glm(final.model, 
                                        newdata = as.data.frame(x.rev),
                                        type = "response"))
  
  # get predicted benefit
  pred.ben.raw <- probs - probs.flipped.w
  pred.ben     <- ifelse(w == 1, -pred.ben.raw, pred.ben.raw)
  
  
  ### 5. housekeeping: make sure the naming is consistent
  # update variable names
  colnames.final.temp <- colnames(x.final)
  colnames.final      <- rep(NA_character_, length(colnames.final.temp))
  
  for(i in 1:p){
    idx <- grepl(pattern = paste0("x", i), x = colnames.final.temp, fixed = TRUE)
    colnames.final[idx] <- gsub(paste0("x", i), colnames.orig[i], colnames.final.temp[idx])
  }
  colnames.final[is.na(colnames.final)] <- "w" # fill up 'w'
  colnames(x.final)                     <- colnames.final
  colnames(x)                           <- colnames.orig
  formula.final.model                   <- paste0("y ~ ", 
                                                  paste(colnames.final, collapse = " + "))
  
  # coefficients
  coefficients <- final.model$coefficients
  names(coefficients) <- c("(Intercept)", colnames.final)
  
  ## 5. fit baseline risk  ----
  # no information on w allowed, so we cannot use the retained variables from the effect modeling
  baseline.mod <- risk.model.stage1(X = x, y = y, alpha = alpha)
  basline.risk <- transform.to.probability(baseline.mod$lp)
  
  ## 6. return ----
  return(list(
    inputs = list(X = x, w = w, y = y),
    baseline.model = baseline.mod,
    effect.model = list(formula = formula.final.model, 
                        selected.data = x.final,
                        glm.obj = final.model,
                        coefficients = coefficients,
                        model.building = list(test.stat = test.stat,
                                              df = df,
                                              pval = pval)),
    risk.regular.w = probs,
    risk.flipped.w = probs.flipped.w,
    risk.baseline = basline.risk,
    predicted.benefit = pred.ben,
    predicted.benefit.raw = pred.ben.raw,
    ate.hat = mean(pred.ben),
    type = "effect.modeling"
  ))
}


# TODO: remove unneccesary functions before execution
