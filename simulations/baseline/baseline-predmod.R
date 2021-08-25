rm(list = ls()) ; cat("\014")

source(paste0(getwd(), "/simulations/baseline/baseline-funs.R"))
source(paste0(getwd(), "/funs/loader.R"))


set.seed(76432)
n <- 1000
p <- 10
theta <- -2.5
time.at.risk.baseline <- 7
survivors.additional.time <- 3

# with binary outcome variable, we can hardly model treatment effect homogeneity!


data <- dgp.baseline(n = n, p = p, 
                     theta = theta,
                     time.at.risk.baseline = time.at.risk.baseline, 
                     survivors.additional.time = survivors.additional.time)

X <- data$X
Y <- data$Y
W <- data$W
lifeyears <- data$time

data$ATE
ARTE0 <- mean(data$p1 / data$p0) # which one is right?
ARTE1 <- mean(data$p1) / mean(data$p0)


# rate ratio
rate.ratio.obj <- rate.ratio(y = Y, w = W, lifeyears = lifeyears) # TODO: is estimate mean(p1)/mean(p0)?

# risk model
risk.model.obj <- risk.modeling(X = X, y = Y, w = W, alpha = 0.5)
risk.model.obj$C.statistics$c.index.outcome.stage1
risk.model.obj$C.statistics$c.index.outcome.stage2
risk.model.obj$models$coefficients.stage2

# effect model
effect.model.obj <- effect.modeling(X = X, y = Y, w = W, alpha = 0.5, interacted.variables = colnames(X))
effect.model.obj$effect.model$summary
effect.model.obj$C.statistics$c.index.outcome
effect.model.obj$C.statistics$c.index.benefit
effect.model.obj$average.treatment.effect

# DML
dml.obj <- dml(X = X, w = W, y = Y)
dml.obj$summary

# GenericML
gml.obj <- GenericML(Z = X, D = W, Y = Y, num.splits = 50, learners.genericML = c("random.forest", "elastic.net"))
gml.obj$VEIN$best.learners$BLP
gml.obj$VEIN$best.learners$GATES

