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
#' @param lifeyears vector of life years. Default is NULL.
#' @param prediction.timeframe vector of the prediction time frame. Default is NULL.
#' @param m number of imputations
#' @param k number of nearest neighbors for donor imputation
#' 
#' @export 
multiple.imputation <- function(X, y, w, 
                                lifeyears = NULL, 
                                prediction.timeframe = NULL, 
                                m = 10, k = 5){
  
  # initialize
  data <- data.frame(y = y, w = w, X = X)
  imputed.datasets <- list()
  
  for(i in 1:m){
    
    # impute data via kNN
    data.imputed <- VIM::kNN(data, k = k, imp_var = FALSE, imp_suffix = FALSE)
    
    imputed.datasets[[paste0("imputed.dataset.", i)]]$X <- data.imputed[,c(3:ncol(data))]
    imputed.datasets[[paste0("imputed.dataset.", i)]]$y <- data.imputed$y
    imputed.datasets[[paste0("imputed.dataset.", i)]]$w <- data.imputed$w
    
  } # FOR
  
  return(imputed.datasets)
} # FUN



# TODO: add decription
# TODO: add imputation routine for lifeyears and prediction.timeframe
grf.modeling.multiple.imputation <- function(X, y, w,
                                             lifeyears = NULL, 
                                             prediction.timeframe = NULL, 
                                             m = 10, k = 5,
                                             num.trees = 2000, ...){
  
  # perform multiple imputation
  imputed.datasets <- multiple.imputation(X = X, y = y, w = w, m = m, k = k)
  
  # initialize
  grf.model.imputed        <- rep(list(NA_real_), m)
  names(grf.model.imputed) <- paste0("imputed.model_", 1:m)
  
  # fit the model m times
  for(i in 1:m){
    
    grf.model.imputed[[i]] <- 
      grf.modeling(X = imputed.datasets[[i]]$X, 
                   y = imputed.datasets[[i]]$y,
                   w = imputed.datasets[[i]]$w, num.trees = num.trees) 
    
  } # FOR m
  
  # return 
  return(grf.model.imputed)
} # FUN


# imp.list is a list of regression outputs
impaccount_regression <- function(imp.list){
  
  # get the relevant statistics
  m     <- length(imp.list)
  T.hat <- rowMeans(sapply(1:m, function(i) imp.list[[i]][,"Estimate"]))
  W.hat <- rowMeans(sapply(1:m, function(i) imp.list[[i]][,"Std. Error"]^2))
  B.hat <- rowSums(sapply(1:m, function(i) (imp.list[[i]][,"Estimate"] - T.hat)^2))/(m-1)
  V.hat <- W.hat + (m+1)/m * B.hat
  
  # inference
  stderr <- sqrt(V.hat)
  z <- T.hat / stderr
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  out <- cbind(T.hat, stderr, z, p)
  rownames(out) <- rownames(imp.list[[1]])
  colnames(out) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  out
  
} # FUN