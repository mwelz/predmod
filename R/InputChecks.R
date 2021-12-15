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