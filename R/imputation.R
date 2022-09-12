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



# C is list of concordance estimates of a given type, 
# SE a list of associated standard errors
impaccount_concordance_single <- function(C, SE)
{
  m <- length(C)
  T_hat <- mean(sapply(1:m, function(i) C[[i]]))
  W_hat <- mean(sapply(1:m, function(i) SE[[i]]^2)) # square since it's a variance
  B_hat <- sum(sapply(1:m, function(i) (C[[i]] - T_hat)^2 )) / (m - 1)
  V_hat <- W_hat + (m+1)/m * B_hat
  
  return(list(estimate = T_hat, stderr = sqrt(V_hat)))
  
} # FUN


# x is a list of concordance objects
#' Account for imputation uncertainty in concordance estimates
#' @param x A list of objects of type \code{"concordance"}
#' @return Imputation accounted concordance estimates
#' 
#' @export
impaccount_concordance <- function(x)
{
  # input check
  impaccount_classcheck(x = x, what = "concordance")
  
  # initialize
  m <- length(x)
  concordance <- list(outcome_baseline = NULL,
                      outcome = NULL,
                      benefit = NULL)

  # check if there is a concordance estimate for baseline risk models
  c_logi  <- all(
    sapply(seq_len(m), function(i) !is.null(x[[i]]$outcome_baseline$estimate )))
  se_logi <- all(
    sapply(seq_len(m), function(i) !is.null(x[[i]]$outcome_baseline$stderr )))
  
  # if yes, then add to list
  if(c_logi && se_logi){
    concordance$outcome_baseline <-
      impaccount_concordance_single(
        C  = lapply(seq_len(m), function(i) x[[i]]$outcome_baseline$estimate),
        SE = lapply(seq_len(m), function(i) x[[i]]$outcome_baseline$stderr)
        )
  } # IF
  
  
  # concordance estimates for C-outcome and C-for-benefit
  concordance$outcome <-
    impaccount_concordance_single(
      C  = lapply(seq_len(m), function(i) x[[i]]$outcome$estimate),
     SE = lapply(seq_len(m), function(i) x[[i]]$outcome$stderr))
  
  concordance$benefit <-
    impaccount_concordance_single(
      C  = lapply(seq_len(m), function(i) x[[i]]$benefit$estimate),
      SE = lapply(seq_len(m), function(i) x[[i]]$benefit$stderr))
  
  # return
  return(concordance)
  
} # FUN



#' void helper function that checks if all elements in a list \code{x}
#' - are of the same class
#' - are of one of the class listed in \code{what}
#' 
#' @noRd
impaccount_classcheck <- function(x, what)
{
  if(!is.list(x))
  {
    stop("x must be a list")
  } # IF
  
  m <- length(x)
  classes <- sapply(seq_len(m), function(i) class(x[[i]]))
  
  if(!identical(length(unique(classes)), 1L))
  {
    stop("All elements in x must be of the same class")
  }
  
  is_inherited <- sapply(seq_len(m), function(i) inherits(x[[i]], what = what))
  
  if(!all(is_inherited))
  {
    stop(paste("The elements in x are not of class:",
               what))
  } # IF
  
} # FUN



#' Account for imputation uncertainty in predmod models
#' 
#' @param x a list of predmod model objects (imputed)
#' @return A list of imputation-accounted stuff
#' 
#' @export
impaccount <- function(x)
{
   # predmod classes
   classes <-  c("risk_model_crss", 
                "risk_model_surv", 
                "effect_model_crss",
                "effect_model_surv")

  # input check
  impaccount_classcheck(x = x, what = classes)
  
  # run correct function
  if(inherits(x[[1L]], what = classes[c(1L,2L)]))
  {
    # case 1: risk model
    impaccount_risk_model(x = x)
  } else{
    # case 2: effect model
    impaccount_effect_model(x = x)
  } # IF
  
} # FUN
