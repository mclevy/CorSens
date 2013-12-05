#' @title Wrapper for L86 (MVN Correlation Mapping)
#'
#' @description This function is a wrapper for the \code{L86} function that formats and passes required information to\code{L86}, which calculates arbitrary-to-MVN distribution correlation coefficient mappings. This function produces a MVN correlation matrix from the correlation of marginal arbitrary distributions (data matrix columns) from \code{data}, or a specified correlation matrix \code{sigma}.
#' 
#' @details This function calculates a MVN correlation matrix (adjusted arbitrary distributions' correlation coefficients), based on either data correlation from \code{data}, or specified parameters \code{k},\code{sigma}, and \code{mu}, and well as \code{distributions}.The function is used internal to the \code{K12} function or may be used independently.
#'
#' @param data A data matrix with columns taken as arbitrary marginal distributions, from which (Pearson) correlation is calculated and mapped to MVN correlation according to each marginal distribution form (specified by \code{distributions}) using \code{L86}.
#' @param distributions A list of distribution names from the list c("unif", "norm", "lnorm", "gamma"), the length of \code{ncol(data)} or \code{k} (in the case of a simulation).
#' @param k \code{NULL} (default) The number of marginal distributions, specified when no data is given, and when \code{sigma} and \code{mu} are provided (such as in a simulation).
#' @param sigma \code{NULL} (default) The marginal distributions' correlation matrix, specified when no data is given, and when \code{k} and \code{mu} are provided (such as in a simulation).
#' @param mu \code{NULL} (default) The marginal distributions' mean vector (of length \code{k}), specified when no data is given, and when \code{k} and \code{sigma} are provided (such as in a simulation).
#' 
#' @references Liu, Pei-Ling, and Armen Der Kiureghian. "Multivariate Distribution Models with Prescribed Marginals and Covariances." Probabilistic Engineering Mechanics 1, no. 2 (June 1986): 105-112. doi:10.1016/0266-8920(86)90033-0.
#' @references Kucherenko, S., S. Tarantola, and P. Annoni. "Estimation of Global Sensitivity Indices for Models with Dependent Variables." Computer Physics Communications 183, no. 4 (April 2012): 937-946. doi:10.1016/j.cpc.2011.12.020.
#'
#' @return matrix A MVN correlation matrix
#' @export
#' @importFrom MASS fitdistr
#' @import stats
#' 
#' 
CorTransform <- function(data, distributions, k=NULL, sigma = NULL, mu=NULL) { 

  # mu is (optional) desired mean vector
  # k is (optional) number of variables when no data is provided
  
  if (is.null(k)){
    n <- ncol(data)
  } else {
    n <- k
  }
  
  combos <- combn(seq(1:n), 2) # combinatins of variables in pairs of 2
  
  # initialize, cor and cov matrix for transformed vars
  Psi_cor <- matrix(NA, ncol=n, nrow=n)
  diag(Psi_cor) <-1
  
  for (i in 1:ncol(combos)){
    xj <- combos[1,i] # 1st pair variable
    xi <- combos[2,i] # 2nd pair variable
    
    if (is.null(mu) == F){
      varj <- mu[xj] # hist(varj)
      vari <- mu[xi]  # hist(vari)
    }else if (is.null(data) == F){
      varj <- data[ ,xj] # hist(varj)
      vari <- data[ ,xi] # hist(vari)
    }
  
    # only used for group 2 distributions (lognormal, gamma)
    delta_j <- sd(varj, na.rm = T)/mean(varj, , na.rm = T) # sigma/mu for distribution j (in Group 2)
    delta_i <- sd(vari, na.rm = T)/mean(vari, , na.rm = T) # sigma/mu for distribution i (in Group 2)
    
    if (is.null(sigma) == T){
      rho_ij <- cor(varj,vari,use="pairwise.complete.obs")
    }else{
      rho_ij <- sigma[xj,xi]
    }
    
    Fm <- L86(delta_j, delta_i, rho_ij, distributions) #  matrix of sigma coefficients (F in Liu, psi in Kucherenko et al. 2012)
    
    Fx <- Fm[distributions[xj],distributions[xi]]
    Psi_cor[xj,xi] <- Psi_cor[xi,xj] <- rho_ij * Fx
    
  }
  return(Psi_cor) # transformed correlation matrix for use with (standard) normal transformed variables
}