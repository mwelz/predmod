# load imputation accounters
source(paste0(getwd(), "/funs/imputation/imputation.R"))

#' applies GRF modeling. Note that risk modeling is not possible here
#' 
#' @param X design matrix, can also be a data frame
#' @param y vector of binary responses. 
#' @param w vector of binary treatment assignments
#' @param lifeyears vector of life years. Default is NULL.
#' @param prediction.timeframe vector of the prediction time frame. Default is NULL.
#' @param num.trees number of trees
#' 
#' @export
grf.modeling <- function(X, y, w,
                         lifeyears = NULL, 
                         prediction.timeframe = NULL, 
                         num.trees = 2000, ...){
  
  # truncate y if necessary
  y.orig    <- y
  
  if(!is.null(lifeyears) & !is.null(prediction.timeframe)){
    lifeyears <- ifelse(lifeyears <= prediction.timeframe, lifeyears, prediction.timeframe) 
    y         <- ifelse(lifeyears <= prediction.timeframe, y, 0)
  } # IF
  
  # get causal forest (for predicted benefit)
  cf <- grf::causal_forest(X = X, Y = y, W = w, num.trees = num.trees, ...)
  
  # causal forest's individual treatment effect estimates are predicted absolute benefit
  predicted.absolute.benefit <- as.numeric(cf$predictions)
  
  # variance estimates 
  predicted.absolute.benefit_variance <- predict(cf, estimate.variance = TRUE)$variance.estimates
  
  # baseline risk
  risk.baseline <- as.numeric(cf$Y.hat)
  
  # ATE
  ate.obj <- grf::average_treatment_effect(cf)
  
  # C statistics
  # match cases based on observed benefit
  matched <- MatchIt::matchit(w ~ pb, data = data.frame(w=w, pb=predicted.absolute.benefit))
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
  pred.ben.abs.paired <- (predicted.absolute.benefit[match.control] +
                            predicted.absolute.benefit[match.treated]) / 2
  c.index.benefit.arr <- Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)
  c.index.benefit     <- list(estimate = unname(c.index.benefit.arr["C Index"]),
                              stderr   = unname(c.index.benefit.arr["S.D."]))
  
  # C index
  c.index.outcome.arr <- Hmisc::rcorr.cens(risk.baseline, y)
  c.index.outcome     <- list(estimate = unname(c.index.outcome.arr["C Index"]),
                              stderr   = unname(c.index.outcome.arr["S.D."]))
  
  # return
  return(list(inputs = list(X = X, w = w, y = y.orig, 
                            lifeyears = lifeyears, 
                            prediction.timeframe = prediction.timeframe, 
                            y.prediction.timeframe = y),
              causal.forest.obj = cf,
              average.treatment.effect = list(estimate = unname(ate.obj["estimate"]),
                                              stderr = unname(ate.obj["std.err"])),
              risk = list(risk.baseline = risk.baseline),
              benefits = list(predicted.absolute.benefit = predicted.absolute.benefit,
                              predicted.absolute.benefit_variance = predicted.absolute.benefit_variance),
              C.statistics = list(c.index.outcome = c.index.outcome,
                                  c.index.benefit = c.index.benefit)))
} # FUN


## helper function for imputation uncertainty in GRF estimates of the absolute benefits
## @param grf.ate.objects = lapply(1:m, function(i) grf.model.imputed[[i]]$benefits)
imputation.accounter_grf.benefits <- function(grf.benefits){
  
  # location estimates
  T.hat <- rowMeans(sapply(1:m, function(i) grf.benefits[[i]]$predicted.absolute.benefit))
  
  # within-imputation variance estimates
  W.hat <- rowMeans(sapply(1:m, function(i) grf.benefits[[i]]$predicted.absolute.benefit_variance))
  
  # between-imputation variance estimates
  B.hat <- rowMeans(sapply(1:m, function(i) grf.benefits[[i]]$predicted.absolute.benefit - T.hat)^2) / (m-1)
  
  return(list(predicted.absolute.benefit = T.hat,
              predicted.absolute.benefit_variance = W.hat + (m+1)/m * B.hat))
  
} # FUN


#' TODO: write documentation
#'
#'
grf.modeling_imputation.accounter <- function(grf.model.imputed){
  
  # initialize
  grf.model.imp.adj <- list()
  m <- length(grf.model.imputed)
  
  # ATE
  grf.model.imp.adj$average.treatment.effect <- 
    imputation.accounter_scalar.location.stderr(lapply(1:m, function(i) grf.model.imputed[[i]]$average.treatment.effect))
  
  # risk
  grf.model.imp.adj$risk$risk.baseline <- 
    imputation.accounter_location(lapply(1:m, function(i) grf.model.imputed[[i]]$risk$risk.baseline))
  
  # benefits
  benefits.ls <- imputation.accounter_grf.benefits(lapply(1:m, function(i) grf.model.imputed[[i]]$benefits))
  
  grf.model.imp.adj$benefits$predicted.absolute.benefit <- 
    benefits.ls$predicted.absolute.benefit
  
  grf.model.imp.adj$benefits$predicted.absolute.benefit_variance <- 
    benefits.ls$predicted.absolute.benefit_variance
  
  # C statistics
  grf.model.imp.adj$C.statistics$c.index.outcome <-
    imputation.accounter_scalar.location.stderr(lapply(1:m, function(i) grf.model.imputed[[i]]$C.statistics$c.index.outcome))
  
  grf.model.imp.adj$C.statistics$c.index.benefit <- 
    imputation.accounter_scalar.location.stderr(lapply(1:m, function(i) grf.model.imputed[[i]]$C.statistics$c.index.benefit))
  
  return(grf.model.imp.adj)
  
} # FUN