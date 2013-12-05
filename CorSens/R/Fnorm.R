#' @title Normal CDF function for \code{K12}
#'
#' @description Function computes normal CDF with data-respective means (data are scaled but not centered); function takes the MVN transformed data from \code{K12}, and returns correlated, uniformly distributed variables.
#' @details The function operates within the \code{K12} function. It returns matrix of uniform, correlated marginal distributions which are subsequently transformed to arbitrary marginal distributions using the \code{Ginv} function. This function is a component of the copula method within the \code{K12} function.
#'
#' @param data A data matrix of MVN transformed data from the \code{K12} function
#'
#' @return matrix
#' @export
#' @import stats
#' 
#' 
Fnorm <- function(data){
  
  Fx <- matrix(NA, ncol=ncol(data), nrow=nrow(data))
  m <- colMeans(data)
  
  for (i in 1:ncol(Fx)){
    Fx[ ,i] <-  pnorm(data[ ,i], mean=m[i]) # hist(Fx[ ,i])
  }
  
  return(Fx)
}