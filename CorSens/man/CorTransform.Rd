\name{CorTransform}
\alias{CorTransform}
\title{Wrapper for L86 (MVN Correlation Mapping)}
\usage{
CorTransform(data = NULL, distributions, k = NULL, sigma = NULL,
  mu = NULL)
}
\arguments{
  \item{data}{A data matrix with columns taken as arbitrary
  marginal distributions, from which (Pearson) correlation
  is calculated and mapped to MVN correlation according to
  each marginal distribution form (specified by
  \code{distributions}) using \code{L86}.}

  \item{distributions}{A list of distribution names from
  the list c("unif", "norm", "lnorm", "gamma"), the length
  of \code{ncol(data)} or \code{k} (in the case of a
  simulation).}

  \item{k}{\code{NULL} (default) The number of marginal
  distributions, specified when no data is given, and when
  \code{sigma} and \code{mu} are provided (such as in a
  simulation).}

  \item{sigma}{\code{NULL} (default) The marginal
  distributions' covariance matrix, specified when no data
  is given, and when \code{k} and \code{mu} are provided
  (such as in a simulation).}

  \item{mu}{\code{NULL} (default) The marginal
  distributions' mean vector (of length \code{k}),
  specified when no data is given, and when \code{k} and
  \code{sigma} are provided (such as in a simulation).}
}
\value{
matrix A MVN correlation matrix
}
\description{
This function is a wrapper for the \code{L86} function that
formats and passes required information to\code{L86}, which
calculates arbitrary-to-MVN distribution correlation
coefficient mappings. This function produces a MVN
correlation matrix from the correlation of marginal
arbitrary distributions (data matrix columns) from
\code{data}, or a specified correlation matrix
\code{sigma}.
}
\details{
This function calculates a MVN correlation matrix (adjusted
arbitrary distributions' correlation coefficients), based
on either data correlation from \code{data}, or specified
parameters \code{k},\code{sigma}, and \code{mu}, and well
as \code{distributions}.The function is used internal to
the \code{K12} function or may be used independently.
}
\references{
Liu, Pei-Ling, and Armen Der Kiureghian. "Multivariate
Distribution Models with Prescribed Marginals and
Covariances." Probabilistic Engineering Mechanics 1, no. 2
(June 1986): 105-112. doi:10.1016/0266-8920(86)90033-0.

Kucherenko, S., S. Tarantola, and P. Annoni. "Estimation of
Global Sensitivity Indices for Models with Dependent
Variables." Computer Physics Communications 183, no. 4
(April 2012): 937-946. doi:10.1016/j.cpc.2011.12.020.
}

