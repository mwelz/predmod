rm(list = ls()) ; cat("\014")
set.seed(123)
n = 100
p = 5

x = list()

for(i in 1:2)
{
  
  X = matrix(runif(n*p), n, p)
  Y = rbinom(n, 1, 0.5)
  W = rbinom(n, 1, 0.5)
  
  x[[i]] = risk_model(X, Y, W)
  
}

rm(list = setdiff(ls(), "x"))

cutoffs = c(0.25, 0.5, 0.75)
baseline_risk = NULL
benefits_risk = FALSE
time_eval = NULL
significance_level = 0.05

#######################
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
    arr_ls$ob_abs[,,i] <- ben$predicted_benefit$absolute[,c("estimate", "stderr")]
    arr_ls$ob_rel[,,i] <- ben$predicted_benefit$relative[,c("estimate", "stderr")]
    arr_ls$or[,,i] <- ben$odds_ratio[,c("estimate", "stderr")]
    
  } # FOR
  
  
  
} # 


