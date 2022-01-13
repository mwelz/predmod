# x is a 3D-array of regression outputs on the 3rd dimension
impaccount_regression_array <- function(x, 
                                        relative = FALSE, 
                                        significance_level = 0.05){
  
  # since x carries 'stderr' in 2nd dimension, but we need the variance, we need to adjust
  x[,2L,] <- x[,2L,]^2L
  
  # get the relevant statistics
  m     <- dim(x)[3L]
  ave   <- apply(x, c(1L,2L), mean)
  T_hat <- ave[, 1L] # location estimate
  W_hat <- ave[, 2L] # within-variance estimate
  
  # between-variance estimate
  B_hat <- rowSums(sapply(1:m, 
                          function(i) (x[,1L,i] - T_hat)^2))/(m-1)
  
  # pooled variance estimate
  V_hat <- W_hat + (m+1)/m * B_hat
  
  # get confidence interval
  stderr <- sqrt(V_hat)
  z      <- stats::qnorm(1 - significance_level / 2)
  
  if(relative){
    cilo <- T_hat * exp(-z * stderr)
    ciup <- T_hat * exp(z * stderr)
  } else{
    cilo <- T_hat - z * stderr
    ciup <- T_hat + z * stderr
  } # IF
  
  # organize output to be consistent with other benefits functions
  out <- cbind(estimate = T_hat,
               ci_lower = cilo,
               ci_upper = ciup,
               stderr = stderr)
  return(out)
  
} # FUN


get_benefits_imputation <- function(x, 
                                    cutoffs = c(0.25, 0.5, 0.75),
                                    baseline_risk = NULL,
                                    benefits_risk = FALSE,
                                    time_eval = NULL,
                                    significance_level = 0.05)
{
  # number of imputation runs
  m <- length(x)
  
  # initialize array of results
  arr <- array(data = NA_real_, 
               dim = c(length(cutoffs) + 1L, 2L, m),
               dimnames = list(NULL, c("estimate", "stderr"), NULL))
  
  # initialize list to store results in
  arr_ls <- list(pb_abs = arr,
                 pb_rel = arr,
                 ob_abs = arr,
                 ob_rel = arr,
                 or     = arr)
  
  for(i in 1:m)
  {
    # get benefits
    ben <- get_benefits(x = x[[i]], 
                        cutoffs = cutoffs,
                        baseline_risk = baseline_risk,
                        time_eval = time_eval,
                        significance_level = significance_level)
    
    # assign results matrices
    arr_ls$pb_abs[,,i] <- ben$predicted_benefit$absolute[,c("estimate", "stderr")]
    arr_ls$pb_rel[,,i] <- ben$predicted_benefit$relative[,c("estimate", "stderr")]
    arr_ls$ob_abs[,,i] <- ben$observed_benefit$absolute[,c("estimate", "stderr")]
    arr_ls$ob_rel[,,i] <- ben$observed_benefit$relative[,c("estimate", "stderr")]
    arr_ls$or[,,i]     <- ben$odds_ratio[,c("estimate", "stderr")]
    
  } # FOR
  
  # account for imputation uncertainty in absolute estimates
  pb_abs <- impaccount_regression_array(arr_ls$pb_abs,
                                        relative = FALSE, 
                                        significance_level = significance_level)
  ob_abs <- impaccount_regression_array(arr_ls$ob_abs,
                                        relative = FALSE, 
                                        significance_level = significance_level)
  
  # account for imputation uncertainty in relative estimates
  pb_rel <- impaccount_regression_array(arr_ls$pb_rel,
                                        relative = TRUE, 
                                        significance_level = significance_level)
  ob_rel <- impaccount_regression_array(arr_ls$ob_rel,
                                        relative = TRUE, 
                                        significance_level = significance_level)
  or     <- impaccount_regression_array(arr_ls$or,
                                        relative = TRUE, 
                                        significance_level = significance_level)
  
  # organize output
  out <- list(predicted_benefit = list(absolute = pb_abs,
                                       relative = pb_rel),
              observed_benefit = list(absolute = ob_abs,
                                      relative = ob_rel),
              odds_ratio = or,
              significance_level = significance_level)
  return(out)

} # FUN

