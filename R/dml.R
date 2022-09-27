
#' Performs DML estimation by using an interactive regression model
#' 
#' @param X a matrix or data frame of covariates
#' @param w a binary vector of treatment status
#' @param status a binary vector of outcomes
#' @param ml_g a regression machine learner; refers to the nuisance function \code{g0(X) = E[Y|X,W]}. Either  'glm', 'random.forest', or 'tree'. Can alternatively be specified by using the mlr3 framework, for example ml_g = mlr3::lrn("regr.ranger", num.trees = 500, mtry = NULL, min.node.size = NULL, max.depth = NULL) for a classification forest, which is also the default. 
#' @param ml_m a classification machine learner; refers to the nuisance function \code{m0(X) = E[W|X]}. Either  'glm', 'random.forest', or 'tree'. Can alternatively be specified by using the mlr3 framework, for example ml_m = mlr3::lrn("classif.ranger", num.trees = 500, mtry = NULL, min.node.size = NULL, max.depth = NULL) for a classification forest, which is also the default.
#' @param significance.level TODO
#' @param n_folds TODO
#' @param score TODO
#' @param trimming_rule TODO
#' @param trimming_threshold TODO
#' @param dml_procedure TODO
#' @param draw_sample_splitting TODO
#' @param apply_cross_fitting TODO
#' 
#' @return Estimates of the causal parameter and a 'DoubleML' object 
#' 
#' @import mlr3verse
#' 
#' @export
dml <- function(X, w, status, 
                ml_g = "random.forest",
                ml_m = "random.forest",
                significance.level = 0.05,
                n_folds = 5,
                score = "ATE",
                trimming_rule = "truncate",
                trimming_threshold = 1e-12,
                dml_procedure = "dml2",
                draw_sample_splitting = TRUE,
                apply_cross_fitting = TRUE){
  
  # surpress messages from mlr3 package during fitting
  lgr::get_logger("mlr3")$set_threshold("warn")
  
  # prepare the machine learner for ml_g
  if(is.environment(ml_g)){
    ml_g <- ml_g
  } else if(ml_g == "glm"){
    
    ml_g <- mlr3::lrn("classif.cv_glmnet", s = "lambda.1se")
    
  } else if(ml_g == "random.forest"){
    
    ml_g <- mlr3::lrn("classif.ranger", num.trees = 500)
    
  } else if(ml_g == "tree"){
    
    ml_g <- mlr3::lrn("classif.rpart")
    
  } else{
    
    stop("Invalid argument for ml_g. Needs to be either 'glm', 'random.forest', 'tree', or an mlr3 object")
    
  } # END IF
  
  
  # prepare the machine learner for ml_m
  if(is.environment(ml_m)){
    ml_m <- ml_m
  } else if(ml_m == "glm"){
    
    ml_m <- mlr3::lrn("classif.cv_glmnet", s = "lambda.1se")
    
  } else if(ml_m == "random.forest"){
    
    ml_m <- mlr3::lrn("classif.ranger", num.trees = 500)
    
  } else if(ml_m == "tree"){
    
    ml_m <- mlr3::lrn("classif.rpart")
    
  } else{
    
    stop("Invalid argument for ml_m. Needs to be either 'glm', 'random.forest', 'tree', or an mlr3 object")
    
  } # END IF
  
  
  # alternative examples for the machine learners:
  # GLM with cross validation: lrn("regr.cv_glmnet",s = "lambda.min") 
  # GLM with fixed lambda: lrn("regr.glmnet", lambda = sqrt(log(n_vars)/(n_obs)))
  # Random Forest (Regression): lrn("regr.ranger")
  # Random Forest (Classification): lrn("classif.ranger")
  # Tree: lrn("regr.rpart")
  
  # matrix interface to DoubleMLData
  dml.data <- DoubleML::double_ml_data_from_matrix(X = X, y = status, d = w)
  
  # specify MLIRM model (estimates ATE if there is treatment effect heterogeneity)
  # -> generalizes PLM
  dml.obj <- DoubleML::DoubleMLIRM$new(dml.data,
                                       ml_g = ml_g,
                                       ml_m = ml_m,
                                       n_folds = n_folds,
                                       n_rep = 1,
                                       score = score,
                                       trimming_rule = trimming_rule,
                                       trimming_threshold = trimming_threshold,
                                       dml_procedure = dml_procedure,
                                       draw_sample_splitting = draw_sample_splitting,
                                       apply_cross_fitting = apply_cross_fitting)
  
  # fit PLR model
  dml.obj$fit(store_predictions = TRUE)
  
  # get CI
  ci <- dml.obj$confint(level = 1 - significance.level)
  cinam <- paste0(paste0(100 * (1 - significance.level), "CI "), c("lower", "upper"))
  
  # evaluate model
  summary <- c(Estimate = dml.obj$all_coef,
               StdError = dml.obj$all_se, 
               CI.lower = ci[1],
               CI.upper = ci[2],
               t.test   = as.numeric(dml.obj$t_stat), 
               p.value  = as.numeric(dml.obj$pval))
  names(summary) <- 
    c("Estimate", "Std. Error", cinam, "t test", "p value")
  
  return(list(summary = summary,
              significance_level = significance.level,
              model = dml.obj))
  
} # END function
