#' @title Arbitrary Inverse CDF for \code{K12}
#'
#' @description The function transforms a matrix of uniform, correlated marginal distributions to correlated arbitrary marginal distributions, using a list of specified parameters or those calculated with \code{ParamFit} function.
#'
#'@details The function operates within the \code{K12} function, and after using the \code{Fnorm} function on MVN transformed data from \code{K12}. It transforms a matrix of uniform, correlated marginal distributions to arbitrary marginal densities. This function is a component of the copula method within the \code{K12} function; the combined use of the \code{Fnorm} function, followed by the \code{Ginv} function on the \code{K12} MVN transformed data is a Gaussian copula.
#'
#' @param data A data matrix of correlated uniform marginal distributions output from the Fnorm function (within the \code{K12} function).
#' @param distributions A vector of distribution names from the list c("unif", "norm", "lnorm", "gamma"), the length of data columns.
#' @param params A list of arbitrary distribution fit parameters; specified in list form, or calculated using the \code{ParamFit} function.
#'
#' @return matrix
#' @export
#' @import stats
#' 
#' 
Ginv <- function(data, distributions, params){
  
  G <- matrix(NA, ncol=ncol(data), nrow=nrow(data))
  distfunc <- paste0("q", distributions)
  
  for (i in 1:ncol(G)){
    distf <- get(distfunc[i]) # extract 
    
    if (distributions[i] == "unif"){
      G[ ,i] <- data[ ,i] # don't transform, already uniform
    } else {
      G[ ,i] <- distf(data[ ,i], params[[i]][1], params[[i]][2]) # hist(G[ ,i])
    }
  }
  
  return(G)
}