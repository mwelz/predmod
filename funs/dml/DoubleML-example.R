rm(list = ls()) ; cat("\014")
# based on https://docs.doubleml.org/r/stable/articles/DoubleML.html

### 0. data ----
# Load bonus data
df_bonus = DoubleML::fetch_bonus(return_type="data.table")

# Simulate data
set.seed(3141)
n_obs = 500
n_vars = 100
theta = 3
X = matrix(rnorm(n_obs*n_vars), nrow=n_obs, ncol=n_vars)
d = X[,1:3]%*%c(5,5,5) + rnorm(n_obs)
y = theta*d + X[, 1:3]%*%c(5,5,5) + rnorm(n_obs)


### 1. The data-backend DoubleMLData ----

# Specify the data and variables for the causal model
dml_data_bonus = DoubleML::DoubleMLData$new(df_bonus,
                                  y_col = "inuidur1",
                                  d_cols = "tg",
                                  x_cols = c("female", "black", "othrace", "dep1", "dep2",
                                             "q2", "q3", "q4", "q5", "q6", "agelt35", "agegt54",
                                             "durable", "lusd", "husd"))

# matrix interface to DoubleMLData
dml_data_sim = DoubleML::double_ml_data_from_matrix(X=X, y=y, d=d)


### 2. Learners to estimate the nuisance models ----
# To estimate our partially linear regression (PLR) model with the double machine learning algorithm, we first have to specify machine learners to estimate m_0 and g_0. For the bonus data we use a random forest regression model and for our simulated data from a sparse partially linear model we use a Lasso regression model.

# surpress messages from mlr3 package during fitting
lgr::get_logger("mlr3")$set_threshold("warn")

learner = mlr3::lrn("regr.ranger", num.trees=500, mtry=floor(sqrt(n_vars)), max.depth=5, min.node.size=2)
ml_g_bonus = learner$clone()
ml_m_bonus = learner$clone()

learner = mlr3::lrn("regr.glmnet", lambda = sqrt(log(n_vars)/(n_obs)))
ml_g_sim = learner$clone()
ml_m_sim = learner$clone()


### 3. Estimate double/debiased machine learning models ----
# We now initialize DoubleMLPLR objects for our examples using default parameters. The models are estimated by calling the fit() method and we can for example inspect the estimated treatment effect using the summary() method. A more detailed result summary can be obtained via the print() method.
set.seed(3141)
obj_dml_plr_bonus = DoubleML::DoubleMLPLR$new(dml_data_bonus, ml_g=ml_g_bonus, ml_m=ml_m_bonus)
obj_dml_plr_bonus$fit()
print(obj_dml_plr_bonus)
bonus.data.estimate <- obj_dml_plr_bonus$summary()

obj_dml_plr_sim = DoubleML::DoubleMLPLR$new(dml_data_sim, ml_g=ml_g_sim, ml_m=ml_m_sim)
obj_dml_plr_sim$fit()
sim.data.estimate <- obj_dml_plr_sim$summary() # truth is theta = 3
