#' Evaluate the rate ratio between two groups
#' 
#' Evaluates the rate ratio between two groups using the rateratio.test package
#' @param status Numeric vector with a unique code for each failure type. Code 0 denotes a censored observation.
#' @param time vector of failure/censoring times.
#' @param w Binary vector of treatment assignment status. Equal to 1 for treatment group and 0 for control group.
#' @param subgroup A integer-valued vector indicating a subgroup
#' @param ... Additional arguments for rateratio.test()
#' @return a rateratio object containing an estimate of the rate ratio between the two groups, as well as the 95% confidence interval and p-value
#' @references  Fay, M. P. (2010). Two-sided exact tests and matching confidence intervals for discrete data. R
#'Journal, 2(1), 53-58. DOI: https://doi.org/10.32614/RJ-2010-008
#' @export
rate_ratio <- function(status, time, w, subgroup = NULL, ...){
  
  # input check
  stopifnot(is.numeric(status))
  stopifnot(is.numeric(time))
  stopifnot(is.numeric(w))
  if(!all(c(0, 1) %in% status)) warning("y is not a binary vector!")
  if(!is.null(subgroup)) stopifnot(is.numeric(subgroup))

  # if no subgroup is specified, all samples are considered
  if(is.null(subgroup)){
    smpl <- seq_len(status)
  } else{
    smpl <- subgroup
  }
  
  # the subgroups, grouped by treatment status
  smpl.w1 <- intersect(which(w == 1), smpl)
  smpl.w0 <- intersect(which(w == 0), smpl)
  
  # rate ratio object
  rr.obj <- rateratio.test::rateratio.test(
    x = c(sum(status[smpl.w1]), sum(status[smpl.w0])),
    n = c(sum(time[smpl.w1]), sum(time[smpl.w0])),
    ...)
  
  # prepare output
  rr <- c(rr.obj$estimate["Rate Ratio"], as.numeric(rr.obj$conf.int), rr.obj$p.value)
  names(rr) <- c("Estimate", "Lower 95%CI", "Upper 95%CI", "p.value")
  
  return(list(rate_ratio = rr,
              rate.ratio_object = rr.obj))
} # FUN