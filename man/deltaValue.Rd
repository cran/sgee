% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deltaValue.R
\name{deltaValue}
\alias{deltaValue}
\title{deltaValue}
\usage{
deltaValue(a, lambda1, lambda2, gamma)
}
\arguments{
\item{a}{Subvector of Estimating Equations that pertains to a particular
group.}

\item{lambda1}{Tuning parameter pertaining to importance of groups.}

\item{lambda2}{Tuning parameter pertaining to importance of individuals.}

\item{gamma}{Scaling / thresholding value.}
}
\value{
The form of a group update for given lambda1, lambda2, and gamma.
}
\description{
Function to calculate an update under a given configuration.
}
\details{
Internal function used by bisee.
}
\note{
Internal function.
}
\author{
Gregory Vaughan
}
\keyword{internal}
