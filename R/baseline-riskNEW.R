baseline_risk <- function(X,
                          status, 
                          alpha = NULL,
                          ...)
{
  
  # input checks
  InputChecks_NA(list(X, status))
  CheckInputs_X(X)
  InputChecks_equal.length2(X, status)
  InputChecks_Y_binary(status)
  
  # assign variable names if there are none
  if(is.null(colnames(X))) colnames(X) <- paste0("V", 1:ncol(X)) 
  
  if(is.null(alpha)){
    
    ## case 1: no regularization
    
    # fit the model
    model.obj <- stats::glm(status~., 
                            family = stats::binomial(link = "logit"), 
                            data = data.frame(status, X), ...)
    
    # get coefficients and linear predictor
    coefs <- as.matrix(model.obj$coefficients, ncol = 1)
    
    # store coefficient matrix as sparse matrix (for consistency with glmnet)
    coefs <- Matrix::Matrix(coefs, sparse = TRUE)
    
    # get linear predictor 
    lp  <- as.numeric(cbind(1,X) %*% coefs)
    
  } else{
    
    ## case 2: regularization
  
    # fit model
    model.obj     <- glmnet::cv.glmnet(X, status, 
                                       family = "binomial", 
                                       alpha = alpha, ...)
    
    # get coefficients at best lambda
    coefs <- glmnet::coef.glmnet(model.obj, s = "lambda.min")
    
    # get indices of kept variables (account for zero indexing)
    kept.vars     <- coefs@i + 1
    
    # get linear predictor 
    lp  <- as.numeric(cbind(1,X)[,kept.vars,drop = FALSE] %*% coefs[kept.vars])
    
  } # IF
  
  
  # return
  return(structure(list(
    risk = stats::plogis(lp),
    linear_predictor = lp,
    coefficients = coefs,
    model = model.obj
  ), class = "baseline_risk"))
  
} # FUN