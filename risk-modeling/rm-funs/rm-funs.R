
quantile.group <- function(x, cutoffs){
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
  kept.vars     <- coefs.obj@i + 1 # adjust for zero-indexing
  coefs         <- coefs.obj@x
  X.star        <- cbind(intercept = 1, X)[,kept.vars]
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
                               intercept = FALSE,
                               offset = oset) 
  
  ### risk prediction:
  ## Note: values aren't bounded by (0,1) because of lp. For true probabilities, set newoffset=0.
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
              risk.pred.regular = as.numeric(probs.regular),
              risk.pred.flipped.w = as.numeric(probs.flipped.w)
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
  
  return(list(
    inputs = list(X = X, w = w, y = y),
    mod.stage1 = stage1$models.stage1.cv,
    mod.stage2 = stage2$mod.stage2,
    linear.predictor = stage1$lp,
    risk.regular.w = stage2$risk.pred.regular,
    risk.flipped.w = stage2$risk.pred.flipped.w,
    predicted.benefit = pred.ben,
    predicted.benefit.raw = pred.ben.raw
  ))
}


transform.to.probability <- function(x) 1 / (1 + exp(-x))


get.benefits <- function(risk.model.obj, cutoffs){
  y <- risk.model.obj$inputs$y
  w <- risk.model.obj$inputs$w
  
  # group observations by their quantile of predicted baseline risk (lp here)
  quantile.groups <- quantile.group(risk.model.obj$linear.predictor, cutoffs)
  
  # extract predicted benefit
  pred.ben <- risk.model.obj$predicted.benefit
  
  ## calculate observed benefit and predicted benefit for each quantile group
  # initialize
  obs.ben.mat           <- as.data.frame(matrix(NA, ncol(quantile.groups), 4))
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
    group.predicted.benefit = pred.ben.mat
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


# TODO: remove unneccesary functions before execution
