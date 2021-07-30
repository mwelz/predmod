#' calculates C index for benefit
#' 
#' @param y binary response vector
#' @param w binary treatment assignment vector
#' @param predicted.benefit vector of predicted benefits
#' @return point estimate and standard error of C index for benefit
#' 
#' @export
C.index.benefit <- function(y, w, predicted.benefit){
  
  # match cases based on observed benefit
  matched <- suppressWarnings(MatchIt::matchit(w ~ pb, data = data.frame(w = w, pb = predicted.benefit)))
  match.treated <- as.numeric(rownames(matched$match.matrix))
  match.control <- as.numeric(matched$match.matrix[,1])
  
  # remove unpaired observations
  no.pairing <- which(is.na(match.control))
  if(length(no.pairing) > 0){
    match.treated <- match.treated[-no.pairing]
    match.control <- match.control[-no.pairing]
  } # IF
  
  # calculate C for benefit by using predicted risk (with regular w)
  obs.ben             <- y[match.control] - y[match.treated]
  pred.ben.abs.paired <- (predicted.benefit[match.control] +
                            predicted.benefit[match.treated]) / 2
  c.index.benefit.arr <- Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)
  
  # return
  return(list(estimate = unname(c.index.benefit.arr["C Index"]),
              stderr   = unname(c.index.benefit.arr["S.D."])))
} # FUN


#' calculates C index outcome
#' 
#' @param y binary response vector
#' @param risk.prediction vector of risk predictions (i.e. Pr(y = 1))
#' @return point estimate and standard error of C index outcome
#' 
#' @export
C.index.outcome <- function(y, risk.prediction){
  
  hmisc.obj <- Hmisc::rcorr.cens(risk.prediction, y)
  return(list(estimate = unname(hmisc.obj["C Index"]),
              stderr   = unname(hmisc.obj["S.D."])))
} # FUN