# required for the function "quantile.group"
source(paste0(getwd(), "/funs/estimation-funs.R"))


#' Estimates the BLP parameters based on the main sample M. 
#' 
#' @param D a binary vector of treatment status of length _|M|_
#' @param Y a vector of responses of length _|M|_
#' @param propensity.scores a vector of propensity scores of length _|M|_
#' @return BLP coefficients with inference statements
#' 
#' @export
#' 
#' TODO: implement same with HT transformation! 
get.BLP.params.classic <- function(D, Y, propensity.scores){
  
  # prepare weights
  weights <- 1 / (propensity.scores * (1 - propensity.scores))
  
  # prepare covariate matrix
  X <- data.frame(B = proxy.baseline, 
                  S = proxy.cate,
                  beta1 = D - propensity.scores, 
                  beta2 = (D - propensity.scores) * (proxy.cate - mean(proxy.cate))) 
  
  # fit weighted linear regression by OLS
  blp.obj <- lm(Y ~., data = data.frame(Y, X), weights = weights)
  
  # extract coefficients
  coefficients     <- summary(blp.obj)$coefficients
  
  # inference on beta2: test the null that it is 1) = 0, 2) = 1.
  beta2.inference <- matrix(NA_real_, 2, 2)
  rownames(beta2.inference) <- c("H0: beta2 = 0", "H0: beta2 = 1")
  colnames(beta2.inference) <- c("t value", "Pr(>|t|)")
  beta2.inference[1,] <- coefficients["beta2", c("t value", "Pr(>|t|)")]
  beta2.inference[2, "t value"] <- 
    (coefficients["beta2", "Estimate"] - 1) / coefficients["beta2", "Std. Error"]  
  beta2.inference[2, "Pr(>|t|)"] <- 
    2 * pt(beta2.inference[2, "t value"], df = blp.obj$df.residual, lower.tail = FALSE)
  
  return(list(lm.obj = blp.obj, 
              blp.coefficients = blp.obj$coefficients[c("beta1", "beta2")],
              coefficients = coefficients,
              beta2.inference = beta2.inference))
  
} # END FUN


#' Estimates the GATES parameters based on the main sample M. 
#' 
#' @param D a binary vector of treatment status of length _|M|_
#' @param Y a vector of responses of length _|M|_
#' @param propensity.scores a vector of propensity scores of length _|M|_
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @return GATES coefficients 
#' 
#' @export
#' 
#' TODO: implement same with HT transformation! 
get.GATES.params.classic <- function(D, Y, 
                                     propensity.scores, 
                                     group.membership.main.sample){
  
  # make the group membership a binary matrix
  groups <- 1 * group.membership.main.sample
  
  # prepare weights
  weights <- 1 / (propensity.scores * (1 - propensity.scores))
  
  # prepare covariate matrix
  X <- data.frame(B = proxy.baseline, 
                  S = proxy.cate,
                  (D - propensity.scores) * groups)
  colnames(X) <- c(colnames(X)[c(1,2)], paste0("gamma", 1:ncol(groups)))
  
  # fit weighted linear regression by OLS
  gates.obj <- lm(Y ~., data = data.frame(Y, X), weights = weights)
  
  # extract coefficients
  coefficients                 <- summary(gates.obj)$coefficients
  gates.coefficients           <- coefficients[paste0("gamma", 1:ncol(groups)), 1]
  gates.coefficients.quantiles <- colnames(groups)
  names(gates.coefficients.quantiles) <- paste0("gamma", 1:ncol(groups))
  
  return(list(lm.obj = gates.obj, 
              gates.coefficients = gates.coefficients,
              gates.coefficients.quantiles = gates.coefficients.quantiles,
              coefficients = coefficients))

} # END FUN


#' Estimates the CLAN parameters in the main sample
#' 
#' @param Z.clan.main.sample a matrix with _|M|_ rows. Each column represents a variable for which CLAN shall be performed.
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @return The two CLAN parameters ("most" affected and "least" affected) for each variable in Z.clan.main.sample
#' 
#' @export
get.CLAN.parameters <- function(Z.clan.main.sample, group.membership.main.sample){
  
  K <- ncol(group.membership.main.sample)
  delta1 <- sapply(1:ncol(Z.clan.main.sample), function(j){
    mean(Z.clan.main.sample[group.membership.main.sample[, 1], j])} )
  deltaK <- sapply(1:ncol(Z.clan.main.sample), function(j){ 
    mean(Z.clan.main.sample[group.membership.main.sample[, K], j])} )
  names(delta1) <- names(deltaK) <- colnames(Z.clan.main.sample)
  return(list(delta.1_parameters = delta1,
              delta.K_parameters = deltaK))
  
} # END FUN


#' returns the two parameters that are used to find the best ML method
#' 
#' @param BLP.obj an object as returned by get.CLAN.parameters()
#' @param GATES.obj an object as returned by get.BLP.parameters()
#' @param proxy.cate.main.sample Proxy CATE estimators for the main sample
#' @param group.membership.main.sample a logical matrix with _M_ rows that indicate 
#' the group memberships (such a matrix is returned by the function quantile.group())
#' @return lambda and lambda.bar parameters
#' 
#' @export
best.ml.method.parameters <- function(BLP.obj,
                                      GATES.obj, 
                                      proxy.cate.main.sample, 
                                      group.membership.main.sample){
  
  return(list(lambda = as.numeric(BLP.obj$blp.coefficients["beta2"]^2 * var(proxy.cate.main.sample)),
              lambda.bar = as.numeric(colSums(group.membership.main.sample) %*%  GATES.obj$gates.coefficients^2)))
  
} # END FUN
