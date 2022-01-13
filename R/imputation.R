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
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  
  out <- cbind(T.hat, stderr, z, p)
  rownames(out) <- rownames(imp.list[[1]])
  colnames(out) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  out
  
} # FUN


# x is a 3D-array of regression outputs on the 3rd dimension
impaccount_regression_array <- function(x){
  
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
  
  # inference
  stderr <- sqrt(V_hat)
  z      <- T_hat / stderr
  p      <- 2.0 * stats::pnorm(abs(z), lower.tail = FALSE)
  
  #### TODO: continue here!
  # BTW, why is there no p-value estimate in get_benefits?
  # do we add this later with an accessor?
  
  out <- cbind(T.hat, stderr, z, p)
  rownames(out) <- rownames(imp.list[[1]])
  colnames(out) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  out
  
} # FUN