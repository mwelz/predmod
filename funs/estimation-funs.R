
transform.to.probability <- function(x) 1 / (1 + exp(-x))


quantile.group <- function(x, cutoffs = c(0.25, 0.5, 0.75), quantile.nam = TRUE){
  # cutoffs are the quantile cutoffs (like c(0.25, 0.5, 0.75))
  q         <- quantile(x, cutoffs)
  q         <- c(-Inf, q, Inf)
  groups    <- as.character(cut(x, breaks = q, include.lowest = FALSE, right = TRUE, dig.lab = 3))
  group.nam <- unique(groups)
  group.nam <- group.nam[order(
    as.numeric(substr(sub("\\,.*", "", group.nam), 2, stop = 1e8L)), 
    decreasing = FALSE)] # ensure the order is correct
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
  
  if(quantile.nam){
    colnames(group.mat) <- nam
  } else{
    colnames(group.mat) <- group.nam
  }
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
  
  # absolute predicted benefit
  pred.ben.abs.raw <- stage2$risk.pred.regular - stage2$risk.pred.flipped.w
  pred.ben.abs     <- ifelse(w == 1, pred.ben.abs.raw, -pred.ben.abs.raw) 
  
  # relative predicted benefit
  pred.ben.rel.raw <- stage2$risk.pred.regular / stage2$risk.pred.flipped.w
  pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)
  
  # coefficients
  coefs.stage1 <- as.matrix(glmnet::coef.glmnet(stage1$models.stage1.cv, s = "lambda.min"))
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
<<<<<<< HEAD
  pred.ben.abs.paired = (pred.ben.abs[control]+pred.ben.abs[treated])/2 #Pairs are matched by predicted benefit; average over the pair
 
  # calculate C for benefit by using predicted risk (with regular w)
   c.index.benefit = unname(Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)[1])
   c.index.stage1.youtcome <- unname(Hmisc::rcorr.cens(stage1$lp, y)[1])
=======
  pred.ben.abs.paired = pred.ben.abs[control]-pred.ben.abs[treated]
  c.index.benefit = unname(Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)[1])
  
  # calculate C for benefit by using predicted risk (with regular w)
  c.index.old <- unname(Hmisc::rcorr.cens(pred.ben.abs, y)[1])
  c.index.youtcome <- unname(Hmisc::rcorr.cens(stage1$lp, y)[1])
>>>>>>> eb1103aba14ce20b24a8624d9ad41704188af9b1
  
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
    predicted.absolute.benefit = pred.ben.abs,
    predicted.absolute.benefit.raw = pred.ben.abs.raw,
    predicted.relative.benefit = pred.ben.rel,
    predicted.relative.benefit.raw = pred.ben.rel.raw,
    ate.hat = mean(pred.ben.abs),
    c.index.old = c.index.old,
    c.index.benefit = c.index.benefit,
<<<<<<< HEAD
    c.index.stage1.youtcome = c.index.stage1.youtcome
=======
    c.index.youtcome = c.index.youtcome
>>>>>>> eb1103aba14ce20b24a8624d9ad41704188af9b1
  ))
}


absolute.observed.benefit <- function(y, w, significance.level = 0.05){
  # calculates the absolute observed benefit, which corresponds to the difference in mean(y[W=w]), along the associated confidence interval.
  
  ttest  <- t.test(y[w==1], y[w==0], paired = FALSE, var.equal = FALSE) 
  aob    <- unname(ttest$estimate[1] - ttest$estimate[2])
  se.aob <- unname(ttest$stderr) 
  
  # quantile of the standard normal distribution
  z <- qnorm(1 - significance.level/2)
  
  return(c(estimate = aob, 
           ci.lower = aob - z * se.aob,
           ci.upper = aob + z * se.aob,
           stderr   = se.aob))
}


relative.observed.benefit <- function(y, w, significance.level = 0.05){
  # calculates the relative observed benefit, which corresponds to a risk ratio, along the associated confidence interval. 
  
  #' contingency table with the following structure
  #'      Y=1 Y=0
  #'  W=1  a   b
  #'  W=0  c   d
  a <- sum(w == 1 & y == 1)
  b <- sum(w == 1 & y == 0)
  c <- sum(w == 0 & y == 1)
  d <- sum(w == 0 & y == 0)
  
  # estimates for the probabilities
  p1.hat <- a / (a + b) # estimates Pr(Y = 1 | W = 1)
  p0.hat <- c / (c + d) # estimates Pr(Y = 1 | W = 0)
  
  # empirical risk ratio (the relative observed benefit) and its standard error
  rob    <- p1.hat / p0.hat
  se.rob <- sqrt(b / (a * (a + b)) + d / (c * (c + d)))
  
  # quantile of the standard normal distribution
  z <- qnorm(1 - significance.level/2)
  
  return(c(estimate = rob, 
           ci.lower = rob * exp(-z * se.rob),
           ci.upper = rob * exp(z * se.rob),
           stderr   = se.rob))
}


odds.ratio <- function(y, w, significance.level = 0.05){
  # calculates the odds ratio, along the associated confidence interval. 
  
  #' contingency table with the following structure
  #'      Y=1 Y=0
  #'  W=1  a   b
  #'  W=0  c   d
  a <- sum(w == 1 & y == 1)
  b <- sum(w == 1 & y == 0)
  c <- sum(w == 0 & y == 1)
  d <- sum(w == 0 & y == 0)
  
  # estimates for the probabilities
  p1.hat <- a / (a + b) # estimates Pr(Y = 1 | W = 1)
  p0.hat <- c / (c + d) # estimates Pr(Y = 1 | W = 0)
  
  # empirical odds ratio and its standard error
  or <- (p1.hat / (1 - p1.hat) ) / (p0.hat / (1 - p0.hat))
  se.or <- sqrt(1/a + 1/b + 1/c + 1/d)
  
  # quantile of the standard normal distribution
  z <- qnorm(1 - significance.level/2)
  
  return(c(estimate = or, 
           ci.lower = or * exp(-z * se.or),
           ci.upper = or * exp(z * se.or),
           stderr   = se.or))
}


predicted.benefit <- function(predicted.benefits, significance.level = 0.05){
  # calculates the predicted benefit, which is, depending on the input, either the relative predicted benefit or the absolute predicted benefit, along the associated confidence interval.
  
  ttest <- t.test(predicted.benefits)
  pb    <- unname(ttest$estimate)
  se.pb <- unname(ttest$stderr)
  
  # quantile of the standard normal distribution
  z <- qnorm(1 - significance.level/2)
  
  return(c(estimate = pb, 
           ci.lower = pb - z * se.pb,
           ci.upper = pb + z * se.pb,
           stderr   = se.pb))
}


get.benefits <- function(pred.model.obj, 
                         cutoffs = c(0.25, 0.5, 0.75),
                         significance.level = 0.05){
  
  # calculates the observed and predicted relative and absolute benefits as well as the odds ratio along the associated confidence intervals.
  
  # extract outcome and treatment status
  y <- pred.model.obj$inputs$y
  w <- pred.model.obj$inputs$w
  
  # group observations by their quantile of predicted baseline risk
  quantile.groups <- quantile.group(pred.model.obj$risk.baseline, cutoffs)
  
  # get predicted benefit (relative and absolute)
  rel.pred.ben <- pred.model.obj$predicted.relative.benefit
  abs.pred.ben <- pred.model.obj$predicted.absolute.benefit
  
  ## calculate observed benefit and predicted benefit for each quantile group
  # initialize
  rel.obs.ben.mat           <- as.data.frame(matrix(NA_real_, ncol(quantile.groups), 5))
  colnames(rel.obs.ben.mat) <- c("quantile", "estimate", "ci.lower", "ci.upper", "stderr")
  rel.obs.ben.mat$quantile  <- gsub(" .*$", "", colnames(quantile.groups))
  abs.obs.ben.mat  <- rel.obs.ben.mat
  abs.pred.ben.mat <- rel.obs.ben.mat
  rel.pred.ben.mat <- rel.obs.ben.mat
  or.mat           <- rel.obs.ben.mat
  
  for(i in 1:ncol(quantile.groups)){
    group <- quantile.groups[,i]
    
    # absolute observed benefit
    abs.obs.ben.mat[i, 2:5] <- 
      absolute.observed.benefit(y[group], w[group], significance.level = significance.level)
    
    # absolute predicted benefit
    abs.pred.ben.mat[i, 2:5] <- 
      predicted.benefit(abs.pred.ben[group], significance.level = significance.level)
    
    # relative observed benefit
    rel.obs.ben.mat[i, 2:5] <- 
      relative.observed.benefit(y[group], w[group], significance.level = significance.level)
    
    # relative predicted benefit
    rel.pred.ben.mat[i, 2:5] <- 
      predicted.benefit(rel.pred.ben[group], significance.level = significance.level)
    
    # odds ratio
    or.mat[i, 2:5] <- 
      odds.ratio(y[group], w[group], significance.level = significance.level)
    
  } # FOR
  
  return(list(absolute.observed.benefit = abs.obs.ben.mat, 
              absolute.predicted.benefit = abs.pred.ben.mat, 
              relative.observed.benefit = rel.obs.ben.mat, 
              relative.predicted.benefit = rel.pred.ben.mat, 
              odds.ratio = or.mat,
              significance.level = significance.level,
              quantile.cutoff.points = cutoffs,
              group.membership = quantile.groups))
}


get.benefits.grf <- function(grf.model.obj, 
                             cutoffs = c(0.25, 0.5, 0.75),
                             significance.level = 0.05){
  # calculates the observed and predicted absolute benefits along the associated confidence intervals for the GRF (relative effects cannot be computed)
  
  # extract outcome and treatment status
  y <- grf.model.obj$inputs$y
  w <- grf.model.obj$inputs$w
  
  # group observations by their quantile of predicted baseline risk
  quantile.groups <- quantile.group(grf.model.obj$risk.baseline, cutoffs)
  
  ## calculate observed benefit and predicted benefit for each quantile group
  # initialize
  abs.obs.ben.mat           <- as.data.frame(matrix(NA_real_, ncol(quantile.groups), 5))
  colnames(abs.obs.ben.mat) <- c("quantile", "estimate", "ci.lower", "ci.upper", "stderr")
  abs.obs.ben.mat$quantile  <- gsub(" .*$", "", colnames(quantile.groups))
  abs.pred.ben.mat <- abs.obs.ben.mat
  
  for(i in 1:ncol(quantile.groups)){
    group <- quantile.groups[,i]
    
    ## absolute observed benefit
    abs.obs.ben.mat[i, 2:5] <- 
      absolute.observed.benefit(y[group], w[group], significance.level = significance.level)
    
    ## absolute predicted benefit
    ate.group <- grf::average_treatment_effect(grf.model.obj$causal.forest.obj, 
                                               subset = group)
    ate <- unname(ate.group["estimate"])
    se  <- unname(ate.group["std.err"]) 
    
    # quantile of the standard normal distribution
    z <- qnorm(1 - significance.level/2)
    
    abs.pred.ben.mat[i, "estimate"] <- ate
    abs.pred.ben.mat[i, "stderr"]   <- se
    abs.pred.ben.mat[i, "ci.lower"] <- ate - z * se
    abs.pred.ben.mat[i, "ci.upper"] <- ate + z * se
    
  } # FOR
  
  return(list(absolute.observed.benefit = abs.obs.ben.mat, 
              absolute.predicted.benefit = abs.pred.ben.mat, 
              significance.level = significance.level,
              quantile.cutoff.points = cutoffs,
              group.membership = quantile.groups))
}



effect.modeling <- function(X, w, y, 
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
  
  ### 1. stage 1: penalized regression on whole set  ----
  # whole set is x.star. We apply sample splitting as suggested in Wasserman and Roeder (2009)
  mod.pm <- glmnet::cv.glmnet(X.star[set1,], y[set1], family = "binomial", alpha = alpha)
  
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
  fit.1     <- glm(y ~., data = data.frame(y, X.kept)[set2,], family = binomial)
  fit.0     <- glm(y ~., data = data.frame(y, X.star.red)[set2,], family = binomial) # reduced model
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
  baseline.mod <- risk.model.stage1(X = X, y = y, alpha = alpha)
  basline.risk <- transform.to.probability(baseline.mod$lp)
  
<<<<<<< HEAD
=======
  
  
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
  pred.ben.abs.paired = pred.ben.abs[control]-pred.ben.abs[treated]
  c.index.benefit = unname(Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)[1])
  
  
  # calculate C index by using predicted risk (with regular w)
  c.index <- DescTools::Cstat(x = probs, resp = y)
>>>>>>> eb1103aba14ce20b24a8624d9ad41704188af9b1
  
  
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
  c.index.youtcome <- unname(Hmisc::rcorr.cens(probs, y)[1])


  ## 6. return ----
  return(list(
    inputs = list(X = X, w = w, y = y),
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
    risk.baseline = basline.risk,
    predicted.absolute.benefit = pred.ben.abs,
    predicted.absolute.benefit.raw = pred.ben.abs.raw,
    predicted.relative.benefit = pred.ben.rel,
    predicted.relative.benefit.raw = pred.ben.rel.raw,
    ate.hat = mean(pred.ben.abs),
<<<<<<< HEAD
    c.index.youtcome = c.index.youtcome,
=======
    c.index = c.index,
>>>>>>> eb1103aba14ce20b24a8624d9ad41704188af9b1
    c.index.benefit =c.index.benefit
  ))
}




group.static <- function(x, group.bounds){
  # if single value supplied: group is all obs with that particular value. 
  # if intervals are supplied: all obs within these bounds (closed set) are a group. Hence, intervals need to be non-overlapping. TODO: add option of having half-closed intervals. 
  # Note that rows that are all FALSE can be produced.
  
  group.mat <- matrix(NA, length(x), length(group.bounds))
  group.nam <- rep(NA_character_, length(group.bounds))
  
  for(i in 1:length(group.bounds)){
    bounds <- group.bounds[[i]]
    
    if(length(bounds) == 1){
      
      # case 1: a single value is supplied
      group.nam[i]   <- as.character(bounds)
      group.mat[, i] <- x == bounds
      
    } else if(length(bounds) == 2){
      
      # case 2: an interval is supplied
      group.nam[i]   <- paste0("[", bounds[1], ",", bounds[2], "]")
      group.mat[, i] <- (bounds[1] <= x) & (x <= bounds[2])
      
    } else stop("Incorrect input. Each element of 'group.bounds' needs to be a scalar or a 2-vector")
  }
  colnames(group.mat) <- group.nam
  return(group.mat)
}


subgroup.plot <- function(pred.model.obj, x, 
                          relative = FALSE,
                          group.bounds = NULL, 
                          quantile.bounds = c(0.25, 0.5, 0.75), 
                          quantile.nam = TRUE,
                          risk.quantile.bounds = c(0.25, 0.5, 0.75)){
  
  # group.bounds <- list(-3.3, c(-3, -2), c(-2+0.0001, 2), c(2 + 0.0001, 3))
  if(!is.null(group.bounds) & !is.null(quantile.bounds)) stop("Either group or quantile!")
  if(!is.null(group.bounds) & !is.list(group.bounds)) stop("group.bounds needs to be a list")
  
  if(!is.null(group.bounds)){
    x.group.mat <- group.static(x, group.bounds = group.bounds)
  } else if(!is.null(quantile.bounds)){
    x.group.mat <- quantile.group(x, cutoffs = quantile.bounds, quantile.nam = quantile.nam)
  } else stop("Please specify quantile bounds or group bounds")
  
  
  # x-axis: risk group
  risk.group.mat <- quantile.group(pred.model.obj$risk.baseline, risk.quantile.bounds) 
  risk.group <- rep(NA_character_, length(x))
  for(nam in colnames(risk.group.mat)){
    risk.group[which(risk.group.mat[,nam])] <- nam
  }
  
  # drop the word 'quantile' and make sure x-axis is in logical order
  risk.group <- factor(gsub( " .*$", "", risk.group))
  lv <- levels(risk.group)
  lv <- lv[order.intervals(lv, quantile.nam = TRUE)]
  risk.group <- factor(risk.group, levels = lv)
  
  ## group by x's values
  x.group <- rep(NA_character_, length(x))
  for(nam in colnames(x.group.mat)){
    x.group[which(x.group.mat[,nam])] <- nam
  }
  
  if(quantile.nam){
    # drop the word 'quantile' and make sure groups are in logical order
    x.group <- factor(gsub( " .*$", "", x.group))
    lv <- levels(x.group)
    lv <- lv[order.intervals(lv, quantile.nam = TRUE)]
    x.group <- factor(x.group, levels = lv)
    leg.tit <- paste0("Quantile of ", deparse(substitute(x)))
  } else{
    x.group <- factor(gsub( " .*$", "", x.group))
    lv <- levels(x.group)
    lv <- lv[order.intervals(lv, quantile.nam = FALSE)]
    x.group <- factor(x.group, levels = lv)
    leg.tit <- paste0("Group of ", deparse(substitute(x)))
  }
  
  # y-axis: predicted benefit
  if(relative){
    pred.ben <- pred.model.obj$predicted.relative.benefit
  } else{
    pred.ben <- pred.model.obj$predicted.absolute.benefit
  }
  
  # prepare data frame
  df <- na.omit(data.frame(pred.ben, risk.group, x.group))
  
  library(ggplot2)
  
  ggplot(data = df, aes(x = risk.group, y = pred.ben, fill = x.group)) + 
    geom_boxplot() +
    theme_bw() +
    xlab("Baseline risk quantile") +
    ylab(ifelse(relative, "Predicted relative benefit", "Predicted absolute benefit")) +
    scale_fill_discrete(name = leg.tit) +
    theme(legend.position = "bottom")
}

order.intervals <- function(intervals, quantile.nam){
  if(quantile.nam){
    is.first <- startsWith(intervals, "<")
    is.last  <- startsWith(intervals, ">")
    cut      <- as.numeric(gsub("[^0-9.-]", "",  gsub(",.*", "", intervals)))
    cut[is.first] <- -Inf
    cut[is.last]  <- Inf
    ord           <- order(cut, decreasing = FALSE)
  } else{
    ord <- order(as.numeric(substr(gsub(",.*", "", intervals), 2, 1e8)), 
                 decreasing = FALSE)
  }
  return(ord)
}



#' Returns a rateratio object as in the package rateratio.test as well as an estimate of the rate ratio.
#' 
#' @param y A binary vector of outcomes
#' @param w A binary vector of treatment status (1 = treatment group)
#' @param lifeyars A vector of life-years
#' @param subgroup A logical vector indicating a subgroup
#' @param ... Additional arguments for rateratio.test()
#' @return a rateratio object and an etsimate of the rate ratio
#' 
#' @export
rate.ratio <- function(y, w, lifeyears, subgroup = NULL, ...){
  
  # input check
  if(!all(c(0, 1) %in% y)) warning("y is not a binary vector!")
  if(any(lifeyears < 0)) warning("Some life years are negative")
  
  # if no subgroup is specified, all samples are considered
  if(is.null(subgroup)){
    smpl <- 1:length(y)
  } else{
    smpl <- subgroup
  }
  
  # the subgroups, grouped by treatment status
  smpl.w1 <- w == 1 & smpl
  smpl.w0 <- w == 0 & smpl
  
  # rate ratio object
  rr.obj <- rateratio.test::rateratio.test(
    x = c(sum(y[smpl.w1]), sum(y[smpl.w0])),
    n = c(sum(lifeyears[smpl.w1]), sum(lifeyears[smpl.w0])),
    ...)
  
  return(list(rate.ratio = unname(rr.obj$estimate["Rate Ratio"]),
              rate.ratio.obj = rr.obj))
}


#' Returns the C index of given risk predictions
#' 
#' @param y A binary vector of outcomes
#' @param risk.predictions A vector of risk predictions for each observation
#' @return The C index
#' 
#' @export
c.index <- function(y, risk.predictions){
  DescTools::Cstat(x = risk.predictions, resp = y)
}



grf.modeling <- function(X, w, y, num.trees = 2000, ...){
  # no relative risk modeling possible!
  # get causal forest (for predicted benefit)
  cf <- grf::causal_forest(X = X, Y = y, W = w, num.trees = num.trees, ...)
  
  # initialize object
  grf.model.obj <- list()
  
  # collect relevant information
  grf.model.obj$inputs            <- list(X = X, w = w, y = y)
  grf.model.obj$causal.forest.obj <- cf
  
  # Y.hat is a random forest's precitions of baseline risk
  grf.model.obj$risk.baseline <- cf$Y.hat
  
  # causal forest's individual treatment effect estimates are predicted absolute benefit
  grf.model.obj$predicted.absolute.benefit <- as.numeric(cf$predictions)
  
  # ATE
  ate.obj <- grf::average_treatment_effect(cf)
  
  # collect everything relevant and return
  grf.model.obj$ate.hat               <- unname(ate.obj["estimate"])
  grf.model.obj$ate.hat.se            <- unname(ate.obj["std.err"])
  
  # TODO: experimental: C statistic (not sure if this is correct as we are using baseline risk)
  grf.model.obj$c.index.youtcome <-  Hmisc::rcorr.cens(rf.model.obj$risk.baseline, y)[1]
 
  
  #Match cases based on observed benefit
  Treatment.formula<- w~grf.model.obj$predicted.absolute.benefit
  matched <- MatchIt::matchit(Treatment.formula)
  treated <- as.numeric(rownames(matched$match.matrix))
  control <- as.numeric(matched$match.matrix[,1])
  
  #Remove unpaired observations
  no.pairing <- which(is.na(matched$match.matrix))
  treated <- treated[-no.pairing]
  control <- control[-no.pairing]
  
  obs.ben <- y[control]-y[treated]
  pred.ben.abs.paired = (grf.model.obj$predicted.absolute.benefit[control]+grf.model.obj$predicted.absolute.benefit[treated])/2
  grf.model.obj$c.index.benefit = unname(Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)[1])
  
  
  #Match cases based on observed benefit
  Treatment.formula<- w~grf.model.obj$predicted.absolute.benefit
  matched <- MatchIt::matchit(Treatment.formula)
  treated <- as.numeric(rownames(matched$match.matrix))
  control <- as.numeric(matched$match.matrix[,1])
  
  #Remove unpaired observations
  no.pairing <- which(is.na(matched$match.matrix))
  treated <- treated[-no.pairing]
  control <- control[-no.pairing]
  
  obs.ben <- y[control]-y[treated]
  pred.ben.abs.paired = grf.model.obj$predicted.absolute.benefit[control]-grf.model.obj$predicted.absolute.benefit[treated]
  grf.model.obj$c.index.benefit = unname(Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)[1])
  
  return(grf.model.obj)
}


calibration.plot <- function( pred.model.obj,
                              cutoffs = c(0.25, 0.5, 0.75), 
                              relative = FALSE,
                              significance.level = 0.05,
                              title = NULL,
                              xlim = NULL,
                              ylim = NULL,
                              flip.sign.of.absolute.benefit = FALSE){
  
  # get observed and predicted benefit by quantile group
  benefits <- get.benefits(pred.model.obj, 
                           cutoffs = cutoffs, 
                           significance.level = significance.level)
  
  # the plot
  library(ggplot2)
  
  # make sure risk quantile is in correct order
  risk.quantile <- factor(benefits$absolute.observed.benefit$quantile)
  lv <- levels(risk.quantile)
  lv <- lv[order.intervals(lv, quantile.nam = TRUE)]
  risk.quantile <- factor(risk.quantile, levels = lv)
  
  # adjust for relative and absolute benefit
  if(relative){
    
    df <- data.frame(pb.means = benefits$relative.predicted.benefit$estimate,
                     ob.means = benefits$relative.observed.benefit$estimate,
                     ob.means.ci.lo = benefits$relative.observed.benefit$ci.lower,
                     ob.means.ci.up = benefits$relative.observed.benefit$ci.upper,
                     risk.quantile = risk.quantile)
  } else{
    
    if(!flip.sign.of.absolute.benefit){
      
      df <- data.frame(pb.means = benefits$absolute.predicted.benefit$estimate,
                       ob.means = benefits$absolute.observed.benefit$estimate,
                       ob.means.ci.lo = benefits$absolute.observed.benefit$ci.lower,
                       ob.means.ci.up = benefits$absolute.observed.benefit$ci.upper,
                       risk.quantile = risk.quantile)
      
    } else{
      
      df <- data.frame(pb.means = -benefits$absolute.predicted.benefit$estimate,
                       ob.means = -benefits$absolute.observed.benefit$estimate,
                       ob.means.ci.lo = -benefits$absolute.observed.benefit$ci.lower,
                       ob.means.ci.up = -benefits$absolute.observed.benefit$ci.upper,
                       risk.quantile = risk.quantile)
      
    }
  }
  
  
  if(is.null(title)){
    title <- paste0("Calibration plot of ", 
                    ifelse(relative, "relative ", "absolute "), "benefit")
  }
  
  ggplot(mapping = aes(x = pb.means,
                       y = ob.means, color = risk.quantile), data = df) +
    geom_point() +
    geom_errorbar(mapping = aes(ymin = ob.means.ci.lo,
                                ymax = ob.means.ci.up)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = ifelse(relative, "Predicted relative benefit", "Predicted absolute benefit"),
         y = ifelse(relative, "Observed relative benefit", "Observed absolute benefit")) +
    theme_light() +
    ggtitle(title) +
    theme(legend.position = "bottom") 
  
}


calibration.plot.grf <- function(grf.model.obj,
                                 cutoffs = c(0.25, 0.5, 0.75), 
                                 significance.level = 0.05,
                                 title = NULL,
                                 xlim = NULL,
                                 ylim = NULL,
                                 flip.sign.of.absolute.benefit = FALSE){
  
  # get observed and predicted benefit by quantile group
  benefits <- get.benefits.grf(grf.model.obj, 
                               cutoffs = cutoffs, 
                               significance.level = significance.level)
  
  # the plot
  library(ggplot2)
  
  # make sure risk quantile is in correct order
  risk.quantile <- factor(benefits$absolute.observed.benefit$quantile)
  lv <- levels(risk.quantile)
  lv <- lv[order.intervals(lv, quantile.nam = TRUE)]
  risk.quantile <- factor(risk.quantile, levels = lv)
  
  
  # retrieve data
  if(!flip.sign.of.absolute.benefit){
    
    df <- data.frame(pb.means = benefits$absolute.predicted.benefit$estimate,
                     ob.means = benefits$absolute.observed.benefit$estimate,
                     ob.means.ci.lo = benefits$absolute.observed.benefit$ci.lower,
                     ob.means.ci.up = benefits$absolute.observed.benefit$ci.upper,
                     risk.quantile = risk.quantile)
    
  } else{
    
    df <- data.frame(pb.means = -benefits$absolute.predicted.benefit$estimate,
                     ob.means = -benefits$absolute.observed.benefit$estimate,
                     ob.means.ci.lo = -benefits$absolute.observed.benefit$ci.lower,
                     ob.means.ci.up = -benefits$absolute.observed.benefit$ci.upper,
                     risk.quantile = risk.quantile)
    
  }
  
  
  if(is.null(title)) title <- "Calibration plot of absolute benefit"
  
  ggplot(mapping = aes(x = pb.means,
                       y = ob.means, color = risk.quantile), data = df) +
    geom_point() +
    geom_errorbar(mapping = aes(ymin = ob.means.ci.lo,
                                ymax = ob.means.ci.up)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = "Predicted absolute benefit",
         y = "Observed absolute benefit") +
    theme_light() +
    ggtitle(title) +
    theme(legend.position = "bottom") 
  
}


# TODO: remove unneccesary functions before execution
