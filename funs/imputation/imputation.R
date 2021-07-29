## helper function to get T.hat for each element in a predictive model object
## 
## @param x = lapply(1:m, function(i) predictive.model.imputed[[i]]$<estimates of interest>)
## @export
imputation.accounter_location <- function(x){
  
  m           <- length(x)
  T.as.matrix <- sapply(1:m, function(i) x[[i]])
  
  if(is.matrix(T.as.matrix)){
    T.hat <- rowMeans(T.as.matrix)
  } else{
    T.hat <- mean(T.as.matrix)
  }
  
  single.obj <- x[[1]]
  
  # make sure return is same object as input
  if(is.matrix(single.obj)){
    
    out <- matrix(T.hat, nrow = nrow(single.obj), ncol = ncol(single.obj), 
                  dimnames = list(rownames(single.obj), colnames(single.obj)))
    
  } else{
    
    out <- T.hat
    names(out) <- names(T.hat)
    
  } # IF
  
  return(out)
  
} # FUN


## helper function for imputation uncertainty scalar estimates with standard errors
## @param x = lapply(1:m, function(i) model.imputed[[i]]$<list of point estimate and standard error>)
imputation.accounter_scalar.location.stderr <- function(x){
  
  m <- length(x)
  
  # mean location estimate
  T.hat <- mean(sapply(1:m, function(i) x[[i]]$estimate))
  
  # average within-imputation variance
  W.hat <- mean(sapply(1:m, function(i) x[[i]]$stderr^2))
  
  # average between-imputation variance
  B.hat <- sum(sapply(1:m, function(i) x[[i]]$estimate - T.hat)^2) / (m-1)
  
  return(
    list(estimate = T.hat, 
         stderr = sqrt(W.hat + (m+1)/m * B.hat))) # pooled variance
} # FUN


#' impute datasets via kNN imputation
#' 
#' @param X matrix of covariates
#' @param y binary outcome variable
#' @param w binary treatment assignment
#' @param m number of imputations
#' @param k number of nearest neighbors for donor imputation
#' 
#' @return 
multiple.imputation <- function(X, y, w, m = 10, k = 5){
  
  # initialize
  data <- data.frame(y = y, w = w, X = X)
  imputed.datasets <- list()
  
  for(i in 1:m){
    
    # impute data via kNN
    data.imputed <- VIM::kNN(data, k = 5, imp_var = FALSE, imp_suffix = FALSE)
    
    imputed.datasets[[paste0("imputed.dataset.", i)]]$X <- data.imputed[,c(3:ncol(data))]
    imputed.datasets[[paste0("imputed.dataset.", i)]]$y <- data.imputed$y
    imputed.datasets[[paste0("imputed.dataset.", i)]]$w <- data.imputed$w
    
  } # FOR
  
  return(imputed.datasets)
} # FUN