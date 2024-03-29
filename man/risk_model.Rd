% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/risk_model_crss.R
\name{risk_model}
\alias{risk_model}
\title{Performs risk modeling by logistic regression.}
\usage{
risk_model(
  X,
  status,
  w,
  z = NULL,
  alpha = 1,
  failcode = 1,
  constant = FALSE,
  LRT = TRUE,
  significance_level = 0.05,
  glm_data = FALSE,
  ...
)
}
\arguments{
\item{X}{Matrix of fixed covariates.}

\item{status}{Numeric vector with a unique code for each failure type. Code 0 denotes a censored observation.}

\item{w}{Binary vector of treatment assignment status. Equal to 1 for treatment group and 0 for control group.}

\item{z}{The \code{z} to be used in the logistic regression model above. If \code{NULL} (default), then the linear predictor of a baseline risk model is used as \code{z}.}

\item{alpha}{The elasticnet mixing parameter for regularization in a potential baseline risk model to obtain the linear predictor to be used as \code{z}. Only applicable if \code{z = NULL}. See \code{\link{baseline_risk}} for details.}

\item{failcode}{Code of status that denotes the failure type of interest. Default is one.}

\item{constant}{Shall \eqn{\beta_3 = 0} be enforced in the logistic regression model above? If \code{TRUE}, then this is equivalent to assuming that there is a constant treatment effect. Default is \code{FALSE}.}

\item{LRT}{shall likelihood ratio test performed to test if there is a constant treatment effect? Default is \code{FALSE}.}

\item{significance_level}{Significance level of possible likelihood ratio test. Only applicable if \code{LRT = TRUE}.}

\item{glm_data}{TODO}

\item{...}{Additional arguments to be passed.}
}
\value{
A \code{predmod_ordinary} object.
}
\description{
The logistic regression model is specified as
\deqn{logit(Pr[y=1|z,w]) = \alpha + \beta_1 w + \beta_2 z + \beta_3 z * w ,}
where \eqn{z} is the typically equal to the linear predictor of a baseline risk model. If \code{constant = TRUE}, it is enforced that \eqn{\beta_3 = 0}.
}
