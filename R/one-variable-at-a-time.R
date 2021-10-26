#' Returns a rateratio object as in the package rateratio.test as well as an estimate of the rate ratio.
#' 
#' @param y A binary vector of outcomes
#' @param w A binary vector of treatment status (1 = treatment group)
#' @param lifeyears A vector of life-years
#' @param prediction.timeframe TODO
#' @param subgroup A logical vector indicating a subgroup
#' @param ... Additional arguments for rateratio.test()
#' @return a rateratio object and an estimate of the rate ratio
#' 
#' @export
rate.ratio <- function(y, w, lifeyears, prediction.timeframe = NULL, subgroup = NULL, ...){
  
  # truncate y if necessary
  if(!is.null(prediction.timeframe)){
    lifeyears <- ifelse(lifeyears <= prediction.timeframe, lifeyears, prediction.timeframe) 
    y         <- ifelse(lifeyears <= prediction.timeframe, y, 0)
  } # IF
  
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
  
  # prepare output
  rr <- c(rr.obj$estimate["Rate Ratio"], as.numeric(rr.obj$conf.int), rr.obj$p.value)
  names(rr) <- c("Estimate", "Lower 95%CI", "Upper 95%CI", "p.value")
  
  return(list(rate.ratio = rr,
              rate.ratio.obj = rr.obj))
} # FUN