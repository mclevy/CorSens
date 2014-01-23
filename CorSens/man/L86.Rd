\name{L86}
\alias{L86}
\title{MVN Correlation Mapping}
\usage{
  L86(delta_j, delta_i, rho_ij, distributions)
}
\arguments{
  \item{delta_j}{A parameter obtained from \code{data} or
  supplied through \code{sigma} and \code{mu} (supplied
  through \code{CorTransform}) equal to the standard
  deviation/mean for marginal distribution 'j' (one of any
  supplied/specified marginal distributions).}

  \item{delta_i}{A parameter obtained from \code{data} or
  supplied through \code{sigma} and \code{mu} (supplied
  through \code{CorTransform}) equal to the standard
  deviation/mean for marginal distribution 'i' (one of any
  supplied/specified marginal distributions).}

  \item{rho_ij}{Correlation between marginal distributions
  'i' and 'j', supplied through \code{CorTransform}.}

  \item{distributions}{A list of distribution names from
  the list c("unif", "norm", "lnorm", "gamma"), supplied
  through \code{CorTransform}.}
}
\value{
  matrix A matrix of mapping parameters "F" (from Liu and
  Kiureghian, 1986) used to generate MVN correlation matrix
  in \code{CorTransform}.
}
\description{
  This function calculates arbitrary-to-MVN distribution
  correlation coefficient mappings, and is used internal to
  the \code{CorTransform} function to generates a MVN
  correlation matrix from the correlation of marginal
  arbitrary distributions (data matrix columns) from
  \code{data}, or a specified correlation matrix
  \code{sigma}.
}
\details{
  This function calculates a MVN correlation matrix
  (adjusted arbitrary distributions' correlation
  coefficients), based on either data covariance from
  \code{data}, or specified parameters
  \code{k},\code{sigma}, and \code{mu}, and well as
  \code{distributions}.The function is used internal to the
  \code{K12} function or may be used independently.
}
\references{
  Liu, Pei-Ling, and Armen Der Kiureghian. "Multivariate
  Distribution Models with Prescribed Marginals and
  Covariances." Probabilistic Engineering Mechanics 1, no.
  2 (June 1986): 105-112. doi:10.1016/0266-8920(86)90033-0.
}

