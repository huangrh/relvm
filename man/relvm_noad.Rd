% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_nonadaptive.R
\name{relvm_noad}
\alias{relvm_noad}
\title{Random Effect Latent Variable Model By Non Adaptive Quadrature Approximation}
\usage{
relvm_noad(object, groups = NULL, fit = list(qpoints = 30, init = NULL,
  predict = TRUE), file = NULL)
}
\arguments{
\item{object}{A mstbl object.}

\item{groups}{A vector of measure group names. The default is NULL, in which
case a vector of all groups will be generated accordingly.}

\item{fit}{A list of fitting parameters. \itemize{ \item qpoints: The numbe
of the quadrature points. \item init: Initial values for mu, fl, and err
term in a list. fl is the factor loading. They will be initialized
generally if it is null. The default is a list with for all mu and one for
others. \item predict: The default is TRUE. \item adaptive: noad, no
adaptive or ad, use adaptive.}}
}
\value{
An list of S3 object of class "relvm" with estimated parametes.
}
\description{
Estimate multiple groups of the random effect latent variable model by gauss
hermite quadrature.
}
\seealso{
\code{\link{mstbl}}
}
