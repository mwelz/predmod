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
grf_model <- function(X, 
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
                        concordance = list(outcome_baseline = C_outcome(y = status_bin, risk = br),
                                           outcome = NULL,
                                           benefit = C_benefit(y = status_bin, 
                                                               w = w,
                                                               pred_ben = benefits$absolute)),
                        models = cf,
                        inputs = list(status = status, status_bin = status_bin,
                                      w = w, failcode = failcode)
  ), 
  class = "grf_ordinary"))
  
} # FUN




#' GRF modeling
#' 
#' @param X Matrix of fixed covariates.
#' @param status Numeric vector with a unique code for each failure type. Code 0 denotes a censored observation.
#' @param time Vector of failure/censoring times.
#' @param w Binary vector of treatment assignment status. Equal to 1 for treatment group and 0 for control group.
#' @param num_trees number of trees
#' @param failcode Code of status that denotes the failure type of interest. Default is one.
#' @param ... Additional arguments
#' 
#' @export
grf_model_survival <- function(X, 
                               status, 
                               time,
                               w,
                               failcode = 1,
                               num_trees = 2000,
                               ...){
  
  # response needs to be binary, so recode status to be binary with 1 = failure due to cause of interest
  status_bin <- ifelse(status == failcode, 1, 0)
  
  cf <- grf::causal_survival_forest(X = X, 
                                    Y = time,
                                    W = w,
                                    D = status_bin,
                                    num.trees = num_trees)
  
  # causal forest's individual treatment effect estimates are predicted absolute benefit
  benefits <- list(absolute = as.numeric(cf$predictions), relative = NULL)
  
  # TODO: how to predict survival here? Amend output
  
  # return
  return(structure(list(benefits = benefits,
                        risk = list(baseline = NULL,
                                    regular = NULL,
                                    counterfactual = NULL),
                        concordance = list(outcome_baseline = NULL,
                                           outcome = NULL,
                                           benefit = NULL),
                        models = cf,
                        inputs = list(status = status, status_bin = status_bin,
                                      w = w, failcode = failcode)
  ), 
  class = "grf_survival"))
  
} # FUN