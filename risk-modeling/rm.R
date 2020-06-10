rm(list = ls()); cat('\014')

source(paste0(getwd(), "/general-funs/general-funs.R"))
source(paste0(getwd(), "/risk-modeling/rm-funs/rm-funs.R"))
set.seed(1)
n             <- 10000
data          <- dgp.screening(n = n, share.lc = 1, reporting.bias = FALSE) # let everybody get LC 
covariates.df <- data$X 
y             <- data$binary.outcomes$outcomes$y
w             <- data$w
X             <- covariates.df[,1:7]
G             <- as.factor(covariates.df[,8])
encoder       <- sufrep::make_encoder(X = X, G = G, method = "one_hot")
X.enc         <- as.matrix(encoder(X, G))
alpha         <- 1
risk.model    <- risk.modeling(X.enc, w, y, alpha, offset.lp = TRUE)
risk.pred     <- transform.to.probability(risk.model$risk.regular.w) # logistic fun to make it a probability

# visualize
hist(risk.model$predicted.benefit, main = "Predicted Benefit")
benefits <- get.benefits(risk.model, cutoffs = c(0.25, 0.5, 0.75))

calibration.plot(risk.model, quantiles = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) 
  
# TODO: put get.benefit and risk.pred in risk.model function
risk.model$mod.stage2.regular.w$beta
risk.model$mod.stage2.flipped.w$beta

