% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correst.R
\name{z_to_corr}
\alias{z_to_corr}
\title{Transforms from z of gaussian copula to correlation of multivariate distribution.}
\usage{
z_to_corr(z, g_cdf, h_cdf)
}
\arguments{
\item{z}{Any real number. The atanh of the underlying gaussian
copula correlation.}

\item{g_cdf}{One marginal cdf. Element i is the probability of
of being less than or equal to i - 1.}

\item{h_cdf}{The other marginal cdf.Element i is the probability of
of being less than or equal to i - 1.}
}
\description{
Transforms from z of gaussian copula to correlation of multivariate distribution.
}
\author{
David Gerard
}
