#' GRF modeling
#' 
#' Estimates a Generalized Random Forest/Causal Forest as described in Wager&Athey (2018)
#' @param X Matrix of fixed covariates.
#' @param status Numeric vector with a unique code for each failure type. Code 0 denotes a censored observation.
#' @param w Binary vector of treatment assignment status. Equal to 1 for treatment group and 0 for control group.
#' @param num_trees number of trees
#' @param failcode Code of status that denotes the failure type of interest. Default is one.
#' @param ... Additional arguments
#' @return A list with predicted absolute benefits, the baseline risk, the causal forest, failure type and failure status
#' @references   Wager S. &  Athey S. (2018) Estimation and Inference of Heterogeneous Treatment Effects using Random Forests,
#'   Journal of the American Statistical Association, 113:523, 1228-1242, DOI: 10.1080/01621459.2017.1319839 
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
                        models = cf,
                        inputs = list(status = status, status_bin = status_bin,
                                      w = w, failcode = failcode)
  ), 
  class = "grf_crss"))
  
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
  class = "grf_surv"))
  
} # FUN
