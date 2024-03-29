% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/risk_model_predict.R
\name{predict.risk_model}
\alias{predict.risk_model}
\title{Predict method for a \code{predmod_crss} object}
\usage{
\method{predict}{risk_model}(object, neww, newz, ...)
}
\arguments{
\item{object}{A \code{predmod_crss} object.}

\item{neww}{Optional treatment assignment variables to calculate benefits with}

\item{newz}{Optional linear predictors}

\item{...}{Additional parameters to be passed down}
}
\value{
A matrix of risk predictions
}
\description{
Predict method for a \code{predmod_crss} object
}
