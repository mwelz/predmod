#' Causal random forest
#' 
<<<<<<< HEAD
#' Estimates a Generalized Random Forest/Causal Forest as described in Wager&Athey (2018)
=======
#' Estimates a causal random forest as described in Wager & Athey (2018)
>>>>>>> master
#' @param X Matrix of fixed covariates.
#' @param status Numeric vector with a unique code for each failure type. Code 0 denotes a censored observation.
#' @param w Binary vector of treatment assignment status. Equal to 1 for treatment group and 0 for control group.
#' @param num_trees number of trees
#' @param failcode Code of status that denotes the failure type of interest. Default is one.
#' @param ... Additional arguments
#' @return A list with predicted absolute benefits, the baseline risk, the causal forest, failure type and failure status
#' @references   Wager S. &  Athey S. (2018) Estimation and Inference of Heterogeneous Treatment Effects using Random Forests,
#'   Journal of the American Statistical Association, 113:523, 1228-1242, DOI: 10.1080/01621459.2017.1319839 
<<<<<<< HEAD
=======
#'   
>>>>>>> master
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


predict.causal_forest <- function(object, newX = NULL, ...)
{
  # get unexported function 'predict.causal_forest'
  predict_cf <- utils::getFromNamespace("predict.causal_forest ", "grf")
  
  # make predictions
  predict_cf(object = object, newdata = newX, ... = ...)
}
