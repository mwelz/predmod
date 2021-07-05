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
  c.index.benefit     <- unname(Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)[1])
  
  # C index
  c.index.outcome <- unname(Hmisc::rcorr.cens(risk.baseline, y)[1])
  
  # return
  return(list(inputs = list(X = X, w = w, y = y.orig, 
                            lifeyears = lifeyears, 
                            prediction.timeframe = prediction.timeframe, 
                            y.prediction.timeframe = y),
              causal.forest.obj = cf,
              average.treatment.effect = list(estimate = unname(ate.obj["estimate"]),
                                              stderr = unname(ate.obj["std.err"])),
              risk = list(risk.baseline = risk.baseline),
              benefits = list(predicted.absolute.benefit = predicted.absolute.benefit),
              C.statistics = list(c.index.outcome = c.index.outcome,
                                  c.index.benefit = c.index.benefit)))
} # FUN
