% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.R
\name{relvm_single_true}
\alias{relvm_single_true}
\title{Estimation Of The Random Effect Latent Variable Model Parameters}
\usage{
relvm_single_true(group, df, inits, predict, use_wt)
}
\arguments{
\item{group}{A measure group name.}

\item{df}{The standardized measure score and measure weight table (alldf).}

\item{inits}{Initial values for mu, fl, and err term in a list. fl is the
factor loading. They will be initialized generally if it is null. The
default is a list with for all mu and one for others.}

\item{predict}{The default is TRUE.}

\item{use_wt}{A vector of strings, possible contains "step1" and "step2",
corresponding to step 1 (optimize the model parameters) and step 2 (predict
the group score) in the latent variable model, respectively.}
}
\value{
An object of S3 class "relvm" with estimated parametes.
}
\description{
Estimate the random effect latent variable model
}
