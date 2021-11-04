#' calculates C index for benefit
#' 
#' @param y binary response vector
#' @param w binary treatment assignment vector
#' @param pred_ben vector of predicted benefits
#' 
#' @return point estimate and standard error of C index for benefit
#' 
#' @export
C_benefit <- function(y, w, pred_ben){
  
  # match cases based on observed benefit
  matched <- suppressWarnings(MatchIt::matchit(w ~ pb, data = data.frame(w = w, pb = pred_ben)))
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
  pred.ben.abs.paired <- (pred_ben[match.control] +
                            pred_ben[match.treated]) / 2
  c.index.benefit.arr <- Hmisc::rcorr.cens(pred.ben.abs.paired, obs.ben)
  
  # return
  return(list(estimate = unname(c.index.benefit.arr["C Index"]),
              stderr   = unname(c.index.benefit.arr["S.D."])))
} # FUN


#' calculates C index outcome
#' 
#' @param y binary response vector
#' @param risk vector of risk predictions (i.e. Pr(y = 1))
#' @return point estimate and standard error of C index outcome
#' 
#' @export
C_outcome <- function(y, risk){
  
  hmisc.obj <- Hmisc::rcorr.cens(risk, y)
  return(list(estimate = unname(hmisc.obj["C Index"]),
              stderr   = unname(hmisc.obj["S.D."])))
} # FUN