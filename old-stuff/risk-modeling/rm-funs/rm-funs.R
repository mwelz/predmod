
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
  pred.ben.abs     <- ifelse(w == 1, -pred.ben.abs.raw, pred.ben.abs.raw)
  
  # relative predicted benefit
  pred.ben.rel.raw <- stage2$risk.pred.regular / stage2$risk.pred.flipped.w
  pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)
  
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
    predicted.absolute.benefit = pred.ben.abs,
    predicted.absolute.benefit.raw = pred.ben.abs.raw,
    predicted.relative.benefit = pred.ben.rel,
    predicted.relative.benefit.raw = pred.ben.rel.raw,
    ate.hat = mean(pred.ben.abs)
  ))
}


get.benefits <- function(pred.model.obj, 
                         cutoffs = c(0.25, 0.5, 0.75),
                         relative = FALSE){
  y <- pred.model.obj$inputs$y
  w <- pred.model.obj$inputs$w
  
  # group observations by their quantile of predicted baseline risk
  quantile.groups <- quantile.group(pred.model.obj$risk.baseline, cutoffs)
  
  # get predicted benefit (relative or absolute, depends on input)
  if(relative){
    pred.ben <- pred.model.obj$predicted.relative.benefit
  } else{
    pred.ben <- pred.model.obj$predicted.absolute.benefit
  }
  
  ## calculate observed benefit and predicted benefit for each quantile group
  # initialize
  obs.ben.mat           <- as.data.frame(matrix(NA_real_, ncol(quantile.groups), 4))
  colnames(obs.ben.mat) <- c("quantile", "mean", "stderr", "df")
  obs.ben.mat$quantile  <- gsub(" .*$", "", colnames(quantile.groups))
  pred.ben.mat          <- obs.ben.mat
  
  for(i in 1:ncol(quantile.groups)){
    group <- quantile.groups[,i]
    
    
    if(!relative){
      ## observed absolute benefit
      # corresponds to difference in mean(y[group & W=w])
      ttest                    <- t.test(y[group & w == 1], y[group & w == 0])
      obs.ben.mat[i, "mean"]   <- unname(ttest$estimate[1] - ttest$estimate[2])
      obs.ben.mat[i, "stderr"] <- unname(ttest$stderr) 
      obs.ben.mat[i, "df"]     <- unname(ttest$parameter)
    } else{
      ## observed relative benefit
      ttest                    <- t.test(y[group & w == 1], y[group & w == 0])
      obs.ben.mat[i, "mean"]   <- mean(y[group & w == 1]) / mean(y[group & w == 0])
      obs.ben.mat[i, "stderr"] <- unname(ttest$stderr) # TODO: wrong, check how it works!
      obs.ben.mat[i, "df"]     <- sum(group) - 1  
    } # IF
    
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


calibration.plot <- function(pred.model.obj,
                             quantiles = c(0.25, 0.5, 0.75), 
                             relative = FALSE,
                             title = NULL,
                             alpha.significance = 0.05){
  
  # get observed and predicted benefit by quantile group
  benefits <- get.benefits(pred.model.obj = pred.model.obj, 
                           cutoffs = quantiles, 
                           relative = relative)
  
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
  
  # make sure risk quantile is in correct order
  risk.quantile <- factor(benefits$group.observed.benefit$quantile)
  lv <- levels(risk.quantile)
  lv <- lv[order.intervals(lv, quantile.nam = TRUE)]
  risk.quantile <- factor(risk.quantile, levels = lv)
  
  df <- data.frame(pb.means = benefits$group.predicted.benefit$mean,
                   ob.means = benefits$group.observed.benefit$mean,
                   ob.means.ci.up = benefits$group.observed.benefit$mean + whisker.obs.ben,
                   ob.means.ci.lo = benefits$group.observed.benefit$mean - whisker.obs.ben,
                   risk.quantile = risk.quantile)
  
  if(is.null(title)) title <- "Calibration plot"
  
  ggplot(mapping = aes(x = pb.means,
                       y = ob.means, color = risk.quantile), data = df) +
    geom_point() +
    geom_linerange(mapping = aes(ymin = ob.means.ci.lo,
                                 ymax = ob.means.ci.up)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = limits, ylim = limits) +
    labs(x = ifelse(relative, "Predicted relative benefit", "Predicted absolute benefit"),
         y = "Observed benefit") +
    theme_classic() +
    ggtitle(title) +
    theme(legend.position = "bottom")
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
  mod.pm <- glmnet::cv.glmnet(x.star[set1,], y[set1], family = "binomial", alpha = alpha)
  
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
  
  # get observed predicted benefit
  pred.ben.abs.raw <- probs - probs.flipped.w
  pred.ben.abs     <- ifelse(w == 1, -pred.ben.abs.raw, pred.ben.abs.raw)
  
  # get relative predicted benefit
  pred.ben.rel.raw <- probs / probs.flipped.w
  pred.ben.rel     <- ifelse(w == 1, pred.ben.rel.raw, 1 / pred.ben.rel.raw)
  
  
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
  
  # summary
  summry <- summary(final.model)$coefficients
  rownames(summry) <- names(coefficients)
  
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
    ate.hat = mean(pred.ben.abs)
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

# TODO: remove unneccesary functions before execution
