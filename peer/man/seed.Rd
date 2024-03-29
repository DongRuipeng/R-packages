% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seed.R
\name{seed}
\alias{seed}
\title{Sequential estimation with eigen-decompsition (SEED)}
\usage{
seed(X, Y, nrank, ntheta = 100)
}
\arguments{
\item{X}{predictors matrix}

\item{Y}{responses matrix}

\item{nrank}{rank of estimator}

\item{ntheta}{number of thresholding parameters for tuning parameter}
}
\value{
\item{C}{estimator for coefficient}
 \item{U}{estimator for left singular vectors}
 \item{V}{estimator for right singular vectors}
 \item{d}{estimator for singular values}
}
\description{
Sequential estimation with eigen-decompsition (SEED)
}
