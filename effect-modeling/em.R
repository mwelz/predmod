rm(list = ls()) ; cat("\014")

# define the logistic function
logistic <- function(x) 1 / (1 + exp(-x))

### 0. Data generation
set.seed(1)
n <- 10000
p <- 5

# treatment assignment (assume RCT)
w <- rbinom(n, 1, 0.5) 

# true treatment effect
tau <- -2

# covariates
x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))

# coefficients (including an intercept)
theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)

# compute Pr(Y = 1 | X) for each individual (with noise)
eps <-  rnorm(n, mean = 0, sd = 0.5)
pi1 <- logistic(tau + as.numeric(cbind(1, x) %*% theta) + eps)
pi0 <- logistic(as.numeric(cbind(1, x) %*% theta) + eps)
hte <- pi1 - pi0

# true ATE
ate <- mean(pi1) - mean(pi0)

# create binary outcomes
y0 <- rbinom(n, 1, pi0)
y1 <- rbinom(n, 1, pi1)
y  <- ifelse(w == 1, y1, y0) # observed outcome

## predictive modeling approach:

effect.modeling <- function(x, w, y, alpha = alpha){
  # split the sample as suggested in Wasserman and Roeder (2009)
  n    <- nrow(x)
  set1 <- sample(1:n, floor(0.5 * n), replace = FALSE)
  set2 <- setdiff(1:n, set1)
  
  # prepare the variable set as in rekkas2019:
  colnames(x) <- paste0("x", 1:p)
  interaction.terms <- sapply(1:p, function(j) ifelse(w == 1, x[,j], 0))
  interaction.terms <- cbind(w, interaction.terms)
  colnames(interaction.terms) <- c("w", paste0("w.", colnames(x)))
  x.star <- cbind(x, interaction.terms)
  alpha <- 1 
  
  # stage 1: penalized regression on whole set x.star
  mod.pm        <- glmnet::cv.glmnet(x.star[set1,], y[set1], family = "binomial", alpha = alpha, 
                                     nfolds = 10)
  
  # stage 2: perform variable selection based on the "best" model
  coefs.obj     <- glmnet::coef.glmnet(mod.pm, s = "lambda.min")
  kept.vars     <- coefs.obj@i 
  if(0 %in% kept.vars){
    kept.vars <- kept.vars[-which(kept.vars == 0)]
  }
  
  
  # stage 3: perform chi-squared test on model with kept variables and other set
  # big model
  x.kept <- x.star[, kept.vars]
  
  # reduced model: (X_lambda, w), where X_lambda are the retained (by the lasso) initial variables 
  x.star.red <- x.star[, 1:(p + 1)]
  x.star.red <- x.star.red[, (kept.vars <= (p + 1))[1:(p+1)]]
  
  # fit both large and reduced model
  fit.1     <- glm(y ~., data = data.frame(y, x.kept)[set2,], family = binomial)
  fit.0     <- glm(y ~., data = data.frame(y, x.star.red)[set2,], family = binomial) # reduced model
  formula.0 <- paste0("y ~ ", paste(colnames(x.star.red), collapse = " + "))
  formula.1 <- paste0("y ~ ", paste(colnames(x.kept), collapse = " + "))
  
  # perform likelihood ratio test
  test.stat <- fit.0$deviance - fit.1$deviance
  df        <- fit.0$df.residual - fit.1$df.residual
  pval      <- pchisq(test.stat, df, lower.tail = FALSE)
  
  # if we reject the null, go with the smaller (the reduced model) to make risk predictions
  if(pval < 0.05){
    final.model         <- fit.0
    x.final             <- x.star.red
    formula.final.model <- formula.0
  } else{
    final.model         <- fit.1
    x.final             <- x.kept
    formula.final.model <- formula.1
  }
  
  # if we reject the null, go with the smaller (the reduced model) to make risk predictions
  probs <- unname(predict.glm(final.model, newdata = as.data.frame(x.final), type = "response"))
  
  return(list(
    formula.final.model = formula.final.model,
    x.final = x.final,
    final.model = final.model,
    probs = probs,
    test.stat = test.stat,
    df = df,
    pval = pval
  ))
}


foo <- effect.modeling(x, w, y, alpha = 1)
  
# TODO: bugfix in risk modeling functions: the kept.vars index at 0 refers to the coefficient, so we do not have zero indexing! (like the way we did it here)
# standardization of X before Lasso; see ESL p. 63
  