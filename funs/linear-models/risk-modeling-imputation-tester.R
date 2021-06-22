rm(list = ls()) ; cat("\014")

# load functions
source(paste0(getwd(), "/funs/linear-models/risk-modeling.R"))
source(paste0(getwd(), "/funs/plotmakers/plotmakers.R"))

#' helper function to get T.hat for each element in a predictive model object
#' 
#' @param x = lapply(1:m, function(i) predictive.model.imputed[[i]]$<estimates of interest>)
#' @export
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


#' TODO: write documentation
#'
#'
imputation.accounter_risk.modeling <- function(predictive.model.imputed){
  
  # initialize
  pred.model.imp.adj <- list()
  
  # stage 1 coefficients
  pred.model.imp.adj$models$coefficients.stage1 <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$models$coefficients.stage1))
  
  # stage 2 coefficients
  pred.model.imp.adj$models$coefficients.stage2 <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$models$coefficients.stage2))
  
  # ATE
  pred.model.imp.adj$average.treatment.effect <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$average.treatment.effect))
  
  # risk regular w
  pred.model.imp.adj$risk$risk.regular.w <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$risk$risk.regular.w))
  
  # risk flipped w
  pred.model.imp.adj$risk$risk.flipped.w <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$risk$risk.flipped.w))
  
  # risk baseline
  pred.model.imp.adj$risk$risk.baseline <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$risk$risk.baseline))
  
  # predicted absolute benefit
  pred.model.imp.adj$benefits$predicted.absolute.benefit <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$benefits$predicted.absolute.benefit))
  
  # predicted relative benefit
  pred.model.imp.adj$benefits$predicted.relative.benefit <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$benefits$predicted.relative.benefit))
  
  # predicted absolute benefit raw
  pred.model.imp.adj$benefits$predicted.absolute.benefit.raw <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$benefits$predicted.absolute.benefit.raw))
  
  # predicted relative benefit raw
  pred.model.imp.adj$benefits$predicted.relative.benefit.raw <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$benefits$predicted.relative.benefit.raw))
  
  # C stat stage 1
  pred.model.imp.adj$C.statistics$c.index.outcome.stage1 <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$C.statistics$c.index.outcome.stage1))
  
  # C stat stage 2
  pred.model.imp.adj$C.statistics$c.index.outcome.stage2 <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$C.statistics$c.index.outcome.stage2))
  
  # C index benefit
  pred.model.imp.adj$C.statistics$c.index.benefit <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$C.statistics$c.index.benefit))
  
  # LP
  pred.model.imp.adj$linear.predictor <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$linear.predictor))
  
  # z
  pred.model.imp.adj$z <- 
    imputation.accounter_location(lapply(1:m, function(i) predictive.model.imputed[[i]]$z))
  
  # return
  return(pred.model.imp.adj)
  
} # FUN


# make data
set.seed(2)
n <- 1000
p <- 5
m <- 10 # number of imputations

predictive.model.imputed <- rep(list(NA_real_), m)
names(predictive.model.imputed) <- paste0("imputed.model_", 1:m)

for(i in 1:m){
  
  # treatment assignment, 50% treatment, 50% control.
  w <- rbinom(n, 1, 0.5) 
  
  # covariates for outcome variable (y)
  X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
  colnames(X) <- paste0("nam", 1:p)
  
  # coefficients (including an intercept)
  theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)
  
  # compute Pr(Y = 1 | X) for each individual (with noise)
  eps <-  rnorm(n, mean = 0, sd = 0.5)
  pi0 <- plogis(as.numeric(cbind(1,X) %*% theta) + eps)
  
  #assume a true relative constant risk reduction of 30%
  scaling <- 0.7
  pi1 <- pi0 * scaling 
  
  # create binary outcomes
  y0 <- rbinom(n, 1, pi0)
  y1 <- rbinom(n, 1, pi1)
  y  <- ifelse(w == 1, y1, y0) # observed outcome
  
  
  # test the risk model
  alpha <- 1
  intercept.stage.2 <- FALSE
  z <- "linear.predictor"
  prediction.timeframe = NULL
  lifeyears = NULL
  
  predictive.model.imputed[[i]] <- 
    risk.modeling(X = X, y = y, w = w, alpha = alpha, 
                  intercept.stage.2 = intercept.stage.2, z = z, 
                  lifeyears = lifeyears, prediction.timeframe = prediction.timeframe)
  
} # FOR imputed datasets


RM.imp <- imputation.accounter_risk.modeling(predictive.model.imputed)



