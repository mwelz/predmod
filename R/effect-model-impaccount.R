#' Account for imputation uncertainty in effect models
#' 
#' @param x a list of effect model objects (imputed)
#' @return A list of imputation-accounted stuff
#' 
#' @export
impaccount_effect_model <- function(x)
{
  
  impaccount_classcheck(x, what = "effect_model_crss")
  
  # number of imputation runs
  m <- length(x)
  
  ## prepare lists
  # benfits
  benefits <- list(absolute = NULL, relative = NULL)
  benefits$absolute <- rowMeans(sapply(1:m, function(i) x[[i]]$benefits$absolute))
  benefits$relative <- rowMeans(sapply(1:m, function(i) x[[i]]$benefits$relative))
  
  # coefficients
  coefficients <- list(baseline = NULL, full = NULL, reduced = NULL)
  
  # account for the case that baseline risk is NULL
  if(all(sapply(1:m, function(i) !is.null(x[[i]]$coefficients$baseline ))))
  {
    
    brimp <- impaccount_baseline_risk(lapply(seq_len(m), function(i) x[[i]]$models$baseline))
    coefficients$baseline <- brimp$coefficients
  }
  
  coefficients$full <- 
    Matrix::Matrix(rowMeans(sapply(1:m, 
                                   function(i) as.matrix(x[[i]]$coefficients$full))),
                   sparse = TRUE)
  
  # jointly retained variables
  var_names <- Reduce(intersect, lapply(1:m, function(i) rownames(x[[i]]$coefficients$reduced)))
  
  arr <- array(NA_real_, dim = c(length(var_names), 4L, m))
  for(i in 1:m) arr[,,i] <- x[[i]]$coefficients$reduced[var_names,]
  tmp <- impaccount_regression_array(x = arr, relative = FALSE)
  z <- tmp[, "estimate"] / tmp[, "stderr"]
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  cf <- cbind(tmp[, "estimate"], tmp[, "stderr"], z, p)
  colnames(cf) <- colnames(x[[1]]$coefficients$reduced)
  rownames(cf) <- var_names
  coefficients$reduced <- cf
  
  # risk
  risk <- list(baseline = NULL, regular = NULL, counterfactual = NULL)
  risk$regular        <- as.matrix(rowMeans(sapply(1:m, function(i) x[[i]]$risk$regular)))
  risk$counterfactual <- as.matrix(rowMeans(sapply(1:m, function(i) x[[i]]$risk$counterfactual)))
  
  if(all(sapply(1:m, function(i) !is.null(x[[i]]$risk$baseline ))))
  {
    risk$baseline     <- as.matrix(rowMeans(sapply(1:m, function(i) x[[i]]$risk$baseline)))
  }
  
  return(list(benefits = benefits,
              coefficients = coefficients,
              risk = risk)) 
} # FOR
