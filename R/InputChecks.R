CheckInputs_X <- function(X) stopifnot(is.matrix(X) & is.numeric(X))


InputChecks_equal.length3 <- function(x, y, z){
  
  X <- as.matrix(x)
  Y <- as.matrix(y)
  Z <- as.matrix(z)
  
  if(!(nrow(X) == nrow(Y) & nrow(X) == nrow(Z) & nrow(Z) == nrow(Y))){
    stop(paste0(deparse(substitute(x)), ", ",
                deparse(substitute(y)), ", ",
                deparse(substitute(z)),
                " need to have an equal number of observations"), call. = FALSE)
  }
} # FUN



InputChecks_equal.length2 <- function(x, y){
  
  X <- as.matrix(x)
  Y <- as.matrix(y)
  
  if(!(nrow(X) == nrow(Y))){
    stop(paste0(deparse(substitute(x)), ", ",
                deparse(substitute(y)),
                " need to have an equal number of observations"), call. = FALSE)
  }
} # FUN



InputChecks_W <- function(W){
  
  # input checks
  if(!(is.numeric(W) & is.vector(W))) stop("W must be a numeric vector", call. = FALSE)
  if((!all(c(0, 1) %in% unique(W))) | (length(unique(W)) != 2)) stop("Treatment assignment W is not binary", call. = FALSE)

} # FUN


InputChecks_Y_binary <- function(Y){
  
  # input checks
  if(!(is.numeric(Y) & is.vector(Y))) stop("Y must be a numeric vector", call. = FALSE)
  if((!all(c(0, 1) %in% unique(Y))) | (length(unique(Y)) != 2)) stop("Treatment assignment Y is not binary", call. = FALSE)
  
} # FUN



InputChecks_NA <- function(list){
  
  if(any(is.na(unlist(list)))) stop("Data contains missing values")
}

InputChecks_newX <- function(newX)
{
  if(!(is.numeric(newX) & is.matrix(newX))) stop("newX must be a numeric matrix", call. = FALSE)
} # FUN


InputChecks_newX_X <- function(newX, object, survival = FALSE)
{
  p <- object$coefficients@Dim[1L] 
  
  if(!survival){
    p <- p - 1L # account for intercept in cross-sectional models
  }
  
  if(!identical(p, ncol(newX))){
    stop(sprintf("The number of variables in newX must be %i", p), 
         call. = FALSE)
  } # IF
  
} # FUN


check_and_adjust_newX <- function(newX, object)
{
  # object must be of class baseline_risk_crss
  nam <- object$covariates
  p <- length(nam)
  
  if(!identical(p, ncol(newX))){
    stop(sprintf("The number of variables in newX must be %i", p), 
         call. = FALSE)
  } # IF
  
  nam_newX <- colnames(newX)
  
  # don't check for column names if newX doesn't have any
  if(!is.null(nam_newX))
  {
    samenames <- nam_newX == nam
    if(!all(samenames))
    {
      stop(paste0("\nInvalid columns detected, namely ", nam_newX[!samenames]), 
           call. = FALSE)
    }
  } else{
    colnames(newX) <- nam
  } # IF
  
  return(newX)
  
} # FUN