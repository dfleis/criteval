% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_cov.R
\name{generate_cov}
\alias{generate_cov}
\title{Generate a covariance matrix with a desired condition number ratio relative to a
correlation matrix}
\usage{
generate_cov(
  kappa.ratio,
  R,
  q = 1,
  .tol.posdef = sqrt(.Machine$double.eps),
  .tol.optimize = sqrt(.Machine$double.eps)
)
}
\arguments{
\item{kappa.ratio}{The target condition number ratio, must be \code{kappa.ratio >= 1}.}

\item{R}{A \eqn{K\times K} strictly positive definite correlation matrix.}

\item{q}{The number of "strong" (high relative variance) features. These features have
their variance set to 1. The remaining \code{K - q} features have variance \code{s < 1}, where
\code{s} is solved for numerically in order to satisfy the condition number ratio. Must
be an integer \verb{1 <= q < K}.}

\item{.tol.posdef}{Tolerance for checking whether the correlation matrix \code{R} is positive definite.}

\item{.tol.optimize}{Tolerance for the \code{optimize} function call.}
}
\value{
A list containing the generated \code{Sigma}, the input \code{R}, the condition
numbers of both, target condition number ratio, and the final observed condition
number ratio given the numerical optimization procedure.
}
\description{
Given a correlation matrix \code{R} and a target condition number ratio \code{kappa.ratio},
this function creates a covariance matrix \code{Sigma} with the property that
\deqn{\frac{\kappa(\Sigma)}{\kappa(R)} = \texttt{kappa.ratio},}
where \eqn{\kappa(A)} denotes the spectral condition number of a matrix \eqn{A}.
}
