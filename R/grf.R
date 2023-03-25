#' GRF modeling
#' 
#' @param X Matrix of fixed covariates.
#' @param status Numeric vector with a unique code for each failure type. Code 0 denotes a censored observation.
#' @param w Binary vector of treatment assignment status. Equal to 1 for treatment group and 0 for control group.
#' @param num_trees number of trees
#' @param failcode Code of status that denotes the failure type of interest. Default is one.
#' @param ... Additional arguments
#' 
#' @export
causal_forest <- function(X, 
                      status, 
                      w,
                      failcode = 1,
                      num_trees = 2000,
                      ...){
  
  # response needs to be binary, so recode status to be binary with 1 = failure due to cause of interest
  status_bin <- ifelse(status == failcode, 1, 0)
  
  # get causal forest (for predicted benefit)
  cf <- grf::causal_forest(X = X, 
                           Y = status_bin,
                           W = w,
                           num.trees = num_trees, ...)
  
  # causal forest's individual treatment effect estimates are predicted absolute benefit
  benefits <- list(absolute = as.numeric(cf$predictions), relative = NULL)
  
  # baseline risk
  br <- as.numeric(cf$Y.hat)
 
  # return
  return(structure(list(benefits = benefits,
                        risk = list(baseline = br,
                                    regular = NULL,
                                    counterfactual = NULL),
                        models = cf,
                        inputs = list(status = status, status_bin = status_bin,
                                      w = w, failcode = failcode)
  ), 
  class = "causal_forest"))
  
} # FUN


#pred <- predict.grf_model(object = x, 
#                          newX = newX)
