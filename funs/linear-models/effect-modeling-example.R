rm(list = ls()) ; cat("\014")

# load functions
source(paste0(getwd(), "/funs/linear-models/effect-modeling.R"))
source(paste0(getwd(), "/funs/plotmakers/plotmakers.R"))


### 0. Data generation ---- 
set.seed(22)
n <- 1000
p <- 5

# treatment assignment, 50% treatment, 50% control.
w <- rbinom(n, 1, 0.5) 

# covariates for outcome variable (y)
X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
colnames(X) <- paste0("nam", 1:p)

# coefficients (including an intercept)
theta <- c(0.2, 0.5, -0.3, 0.7, -0.1, 0.4)

# compute Pr(Y = 1 | X) for each individual (with noise)
pi0 <- plogis(as.numeric(cbind(1,X) %*% theta))
pi1 <- plogis(as.numeric(cbind(1,X) %*% theta) + 1)

# create binary outcomes
y0 <- rbinom(n, 1, pi0)
y1 <- rbinom(n, 1, pi1)
y  <- ifelse(w == 1, y1, y0) # observed outcome


### 1. effect modeling ----
# arguments
alpha = 1
sig.level = 0.05
interacted.variables = colnames(X) 
retained.variables = NULL # A string array of variable names in Z that shall always be retained (also works on interaction variables). If NULL, then no restriction applies. Note that treatment assignment w will always be retained by the function.
significance.level = 0.05
prediction.timeframe = NULL
lifeyears = NULL

# run effect model
em <- effect.modeling(X = X, y = y, w = w, 
                       interacted.variables = interacted.variables, 
                       alpha = alpha, lifeyears = lifeyears, 
                       prediction.timeframe = prediction.timeframe,
                       retained.variables = retained.variables, 
                       significance.level = significance.level)

##############################################################################################

# no information on w allowed, so we cannot use the retained variables from the effect modeling
baseline.mod  <- baseline.risk(X = X, y = y, alpha = 1)
baseline.risk <- baseline.mod$response 


### 2. perform variable selection ----
# We use the strategy in Wasserman and Roeder (2009; Annals of Statistics)
vs.obj <- variable.selection(X = X, y = y, w = w, 
                             interacted.variables = interacted.variables, 
                             alpha = alpha, 
                             retained.variables = retained.variables, 
                             significance.level = significance.level)



# get design matrix for the modeling
if(is.null(interacted.variables)){
  Z <- data.frame(w = w, X)
} else{
  Z <- data.frame(w = w, X,
                  get.interaction.terms.matrix(X = X, 
                                               w = w, 
                                               interacted.variables = interacted.variables))
} # IF

# get Z index of forcefully retained variables
retained.variables.Z.idx <- get.Z.index.of.retained.variables(retained.variables = retained.variables, 
                                                              X = X, Z = Z)

# randomly split data into three roughly equally sized samples
set1 <- sample(1:n, floor(n / 3), replace = FALSE)
set2 <- sample(setdiff(1:n, set1), floor(n / 3), replace = FALSE)
set3 <- setdiff(1:n, c(set1, set2))

# prepare cross-validation: loop over the lambdas of the glmnet path
lambdapath <- get.lambdapath(x = Z[set1,], y = y[set1])

suite <- lapply(1:length(lambdapath), function(...) list(retained.variables = NA, lasso.estimates = NA))
loss  <- rep(NA_real_, length(lambdapath))

for(i in 1:length(lambdapath)){
  
  # penalized logistic regression on first set
  mod <- glmnet::glmnet(Z[set1,], y[set1], 
                        family = "binomial", 
                        alpha = alpha,
                        lambda = lambdapath[i])
  
  # obtain retained variables (excluding intercept)
  kept.vars     <- glmnet::coef.glmnet(mod)@i 
  if(0 %in% kept.vars){
    kept.vars <- kept.vars[-which(kept.vars == 0)]
  } # IF
  
  # make sure that all variables that are forced to be retained will be retained
  kept.vars <- unique(sort(c(kept.vars, retained.variables.Z.idx), decreasing = FALSE))
  
  # add to suite
  suite[[i]]$retained.variables <- kept.vars
  suite[[i]]$lasso.estimates    <- c(as.numeric(mod$a0), as.numeric(mod$beta))
  
  #  calculate logistic least squares estimator with retained variables on first set
  df <- data.frame(y = y[set1], Z[set1, kept.vars])
  colnames(df) <- c("y", colnames(Z)[kept.vars])
  mod <- glm(y~., family =  binomial(link = "logit"), data = df)
  
  # predict Pr(Y=1) on second set
  df <- data.frame(Z[set2, kept.vars])
  colnames(df) <- colnames(Z)[kept.vars]
  p.hat <- predict.glm(mod, newdata = df, type = "response")
  
  # cross-validation: evaluate the loss (cross-entropy here) on the second set
  crossentropy <- rep(NA_real_, length(set2))
  crossentropy[y[set2] == 1] <- -log(p.hat[y[set2] == 1])
  crossentropy[y[set2] == 0] <- -log(1 - p.hat[y[set2] == 0])
  loss[i] <- mean(crossentropy)
  
} # FOR

# find the lambda that minimizes the empirical loss and its associated model
lambda.min <- lambdapath[which.min(loss)]
S.hat      <- suite[[which.min(loss)]]$retained.variables
Z.lambda.min <- Z[, S.hat]

# on third set: use S.hat to calculate logistic least squares estimator with Z.lambda.min
df <- data.frame(y = y[set3], Z.lambda.min[set3,])
colnames(df) <- c("y", colnames(Z.lambda.min))
mod <- glm(y~., family =  binomial(link = "logit"), 
           data = df)

# on third set: hypothesis testing
coeffs <- summary(mod)$coefficients
coeffs <- cbind(coeffs, rep(NA_real_, length(S.hat) + 1))
colnames(coeffs) <- c("Estimate", "Std. Error", "z value", "critical value", "retain?")
critval <- qnorm(significance.level / 2 * length(S.hat), lower.tail = FALSE)
coeffs[, "critical value"] <- critval
coeffs[, "retain?"] <- 1 * (abs(coeffs[, "z value"]) > critval)

# final model
D.hat   <- which(coeffs[-1, "retain?"] == 1)

# make sure that all variables that are forced to be retained will be retained
D.hat <- unique(sort(c(D.hat, retained.variables.Z.idx), decreasing = FALSE))

# get the regularized estimates at the minimizing lambda
regularized.estimates_lambda.min <- suite[[which.min(loss)]]$lasso.estimates
names(regularized.estimates_lambda.min) <- c("(Intercept)", colnames(Z))



# get design matrix associated with the final model
Z     <- vs.obj$final.model$design.matrix_selected.variables

### 3. fit the final effect model and use it for risk estimates ----
# fit the final model, this time on full sample (no penalty required, selection already took place!)
final.model        <- glm(y ~., data = data.frame(y, Z), 
                          family =  binomial(link = "logit"))
predicted.benefits <- effect.model.predicted.benefits(X = X, y = y, w = w, 
                                                      final.model = final.model, Z = Z)




##############################################################################################
calibration.plot(em)
subgroup.plot(em, X[,1])

em$effect.model$summary
em$average.treatment.effect
em$C.statistics
