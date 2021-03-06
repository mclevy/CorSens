#' @title Sensitivity Analysis for Dependent Data
#'
#' @description CorSens is a R package used to do Monte Carlo estimation of sensitivity indices for models with dependent variables (after Kucherenko et al. 2012).
#' 
#' @details When there is no correlation present between model input variables, this method supplies individual and total sensitivity index results roughly equivalent (with differences due to sampling) to the Saltelli sensitivity analysis scheme (2002) (available in the \code{sensitivity} package by Pujol et al.).
#' 
#' The primary function \code{\link{K12}} executes a variance decomposition based sensitivity analysis for models with dependent variables that uses Monte Carlo estimation methods and a Gasussian copula based approach for sampling from arbitrary multivariate distributions. The \code{\link{K12}} funtions allows for flexible input of either data or simulation-based parameters (e.g. \code{mu,sigma}) for model inputs.
#' 
#' Other than the main function \code{\link{K12}}, other functions (used internal to K12) include:
#' 
#' \itemize{
#'  \item \code{\link{ParamFit}} extracts empirical distribution parameter fits (from data).
#'  \item \code{\link{L86}} maps arbitrary-to-MVN distribution correlation coefficients.
#'  \item \code{\link{CorTransform}} is a wrapper for \code{\link{L86}}.
#'  \item \code{\link{Ginv}} calculates arbitrary distribution inverse CDF as a part of the copula function internal to \code{\link{K12}}.
#'  \item \code{\link{Fnorm}} calculates normal CDF as a part of the copula function internal to \code{\link{K12}}.
#'  }
#' 
#' @references Kucherenko, S., S. Tarantola, and P. Annoni. "Estimation of Global Sensitivity Indices for Models with Dependent Variables." Computer Physics Communications 183, no. 4 (April 2012): 937-946. doi:10.1016/j.cpc.2011.12.020.
#' @docType package
#' @name CorSens
NULL