#' calculates C index for benefit
#' 
#' Calculates the C index for benefit as proposed by van Klaveren et al (2018)
#' @param y binary response vector
#' @param w binary treatment assignment vector
#' @param pred_ben vector of predicted benefits
#' 
#' @return point estimate and standard error of C index for benefit
#' 
#' @references van Klaveren D., Steyerberg E.W., Seeruys P.W., Kent D.M. (2018).
#' The proposed concordance-statistic for benefit provided 
#' a useful metric when modeling heterogeneous treatment effects 
#'J Clin Epidemiol. 2018 Feb;94:59-68. doi: 10.1016/j.jclinepi.2017.10.021. Epub 2017 Nov 11. PMID: 29132832; PMCID: PMC7448760.
#' @export
C_benefit <- function(y, w, pred_ben){
  
  InputChecks_equal.length3(y, w, pred_ben)
  
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
  
  loc <- unname(c.index.benefit.arr["C Index"])
  sd <- unname(c.index.benefit.arr["S.D."])
  whisk <- stats::qnorm(0.975) * sd
  
  return(list(estimate = loc,
              stderr   = sd,
              ci95     = c(loc - whisk, loc + whisk)))
} # FUN


#' calculates C index outcome
#' 
#' @param y binary response vector
#' @param risk vector of risk predictions (i.e. Pr(y = 1))
#' @return point estimate and standard error of C index outcome
#' 
#' @export
C_outcome <- function(y, risk){
  
  InputChecks_equal.length2(y, risk)
  
  hmisc.obj <- Hmisc::rcorr.cens(risk, y)
  loc <- unname(hmisc.obj["C Index"])
  sd <-  unname(hmisc.obj["S.D."])
  whisk <- stats::qnorm(0.975) * sd
  return(list(estimate = loc,
              stderr   = sd,
              ci95     = c(loc - whisk, loc + whisk)))
} # FUN



#' Function to predict concordance
#' @param x A predmod object of a risk or effect model
#' @param baseline_risk the baseline risk estimated by the modle
#' @param risk the final risk estimate of the model
#' @param benefit vector of estimated benefits
#' @param status vector of binary outcomes
#' @param w a binary vector of treatment status
#' @export
concordance <- function(x, 
                        baseline_risk = NULL, 
                        risk = NULL, 
                        benefit = NULL,
                        status = NULL,
                        w = NULL ){
  
  stopifnot(inherits(x, what = c("risk_model_crss", 
                                 "risk_model_surv", 
                                 "effect_model_crss",
                                 "effect_model_surv")))
  
  if(is.null(status))
  {
    status0 <- as.numeric(x$inputs$status_bin)
  } else
  {
    status0 <- as.numeric(status)
  }
  
  if(is.null(w))
  {
    w0 <- x$inputs$w
  } else
  {
    w0 <- w
  }
  
  
  ## 1. C baseline
  if(is.null(baseline_risk) && is.null(x$risk$baseline))
  {
    outcome_baseline <- NULL
  } else
  {
    if(is.null(baseline_risk)){
      br <- x$risk$baseline
    } else{
      br <- baseline_risk
    }
    
    # calculate C
    outcome_baseline <- 
      C_outcome(y    = status0, 
                risk = as.numeric(br))
  } # IF
  
  
  ## 2. C outcome
  if(is.null(risk) && is.null(x$risk$regular))
  {
    Coutcome <- NULL
  } else
  {
    if(is.null(risk))
    {
      risk0 <- as.numeric(x$risk$regular)
    } else
    {
      risk0 <- as.numeric(risk)
    } # IF
    
    Coutcome <- C_outcome(y = status0,
                          risk = risk0)
    
  } # IF
  
  
  
  ## 3. predicted absolute benefit for C benefit
  if(is.null(benefit))
  {
    benefit0 <- as.numeric(x$benefits$absolute)
  } else
  {
    benefit0 <- as.numeric(benefit)
  } # IF
  
  Cbenefit <- C_benefit(y = status0,
                        w = w0, 
                        pred_ben = benefit0)
  
 
  return(structure(
    list(
    outcome_baseline = outcome_baseline,
    outcome = Coutcome,
    benefit = Cbenefit
    ), class = "concordance"))
  
} # FUN