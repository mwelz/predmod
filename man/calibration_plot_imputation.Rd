% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-impaccount.R
\name{calibration_plot_imputation}
\alias{calibration_plot_imputation}
\title{makes an imputation-accounted calibration plot for imputed prediction models}
\usage{
calibration_plot_imputation(
  x,
  cutoffs = c(0.25, 0.5, 0.75),
  breaks = NULL,
  relative = FALSE,
  baseline_risk = NULL,
  significance_level = 0.05,
  title = NULL,
  xlim = NULL,
  ylim = NULL,
  flip_sign = FALSE,
  newX = NULL,
  newstatus = NULL,
  neww = NULL,
  newz = NULL,
  shrunk = FALSE
)
}
\arguments{
\item{x}{A list of prediction model objects}

\item{cutoffs}{the cutoff points of quantiles that shall be used for GATES grouping. Default is \code{c(0.25, 0.5, 0.75)}, which corresponds to the quartiles.}

\item{breaks}{List of breaks along which to perform the grouping. If passed, overrules the grouping implied by \code{cutoffs}}

\item{relative}{logical. If \code{TRUE}, then relative benefits will be plotted. Default is \code{FALSE}}

\item{baseline_risk}{A list of baseline risks that shall be used for grouping. If \code{NULL} (default), then the baseline risks as in \code{x} are used.}

\item{significance_level}{significance level for the confidence intervals. Default is 0.05}

\item{title}{optional title of the plot}

\item{xlim}{limits of x-axis}

\item{ylim}{limits of y-xcis}

\item{flip_sign}{logical. Shall the sign of the benefits be flipped?}

\item{newX}{Optional covariate matrix to calculate benefits with}

\item{newstatus}{Optional list of target variables to calculate benefits with}

\item{neww}{Optional list of treatment assignment variables to calculate benefits with}

\item{newz}{Optional linear predictors}

\item{shrunk}{TODO}
}
\description{
makes an imputation-accounted calibration plot for imputed prediction models
}
