#' Account for imputation uncertainty in baseline risk models
#' 
#' @param x a list of baseline risk model objects (imputed)
#' @return A list of imputation-accounted stuff
#' 
#' @export
impaccount_baseline_risk <- function(x) 
{
  m <- length(x)
  risk <- rowMeans(sapply(seq_len(m), function(i) as.numeric(x[[i]]$risk)))
  lp <- rowMeans(sapply(seq_len(m), function(i) as.numeric(x[[i]]$linear_predictor)))
  coefs_full <- rowMeans(sapply(seq_len(m), function(i) as.numeric(x[[i]]$coefficients$full)))
  coefs_full <- Matrix::Matrix(coefs_full, sparse = TRUE)
  dimnames(coefs_full) <- dimnames(x[[1L]]$coefficients$full)
  
  ## coefficient matrix for reduced model
  # jointly retained variables
  var_names <- Reduce(intersect, lapply(seq_len(m), function(i) rownames(x[[i]]$coefficients$reduced)))
  
  arr <- array(NA_real_, dim = c(length(var_names), 4L, m))
  for(i in 1:m) arr[,,i] <- x[[i]]$coefficients$reduced[var_names,]
  tmp <- impaccount_regression_array(x = arr, relative = FALSE)
  z <- tmp[, "estimate"] / tmp[, "stderr"]
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  coefs_red <- cbind(tmp[, "estimate"], tmp[, "stderr"], z, p)
  colnames(coefs_red) <- colnames(x[[1]]$coefficients$reduced)
  rownames(coefs_red) <- var_names
  
  return(
    list(
      risk = as.matrix(risk), linear_predictor = lp, 
      coefficients = list(
        full = coefs_full, reduced = coefs_red
      )
    )
  )
  
} # FUN
 