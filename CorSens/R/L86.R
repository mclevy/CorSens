#' @title MVN Correlation Mapping
#'
#' @description This function calculates arbitrary-to-MVN distribution correlation coefficient mappings, and is used internal to the \code{CorTransform} function to generates a MVN correlation matrix from the correlation of marginal arbitrary distributions (data matrix columns) from \code{data}, or a specified correlation matrix \code{sigma}.
#' 
#' @details This function calculates a MVN correlation matrix (adjusted arbitrary distributions' correlation coefficients), based on either data covariance from \code{data}, or specified parameters \code{k},\code{sigma}, and \code{mu}, and well as \code{distributions}.The function is used internal to the \code{K12} function or may be used independently.
#'
#' @param delta_j A parameter obtained from \code{data} or supplied through \code{sigma} and \code{mu} (supplied through \code{CorTransform}) equal to the standard deviation/mean for marginal distribution 'j' (one of any supplied/specified marginal distributions).
#' @param delta_i A parameter obtained from \code{data} or supplied through \code{sigma} and \code{mu} (supplied through \code{CorTransform}) equal to the standard deviation/mean for marginal distribution 'i' (one of any supplied/specified marginal distributions).
#' @param rho_ij Correlation between marginal distributions 'i' and 'j', supplied through \code{CorTransform}.
#' @param distributions A list of distribution names from the list c("unif", "norm", "lnorm", "gamma"), supplied through \code{CorTransform}.
#' 
#' @references Liu, Pei-Ling, and Armen Der Kiureghian. "Multivariate Distribution Models with Prescribed Marginals and Covariances." Probabilistic Engineering Mechanics 1, no. 2 (June 1986): 105-112. doi:10.1016/0266-8920(86)90033-0.
#'
#' @return matrix A matrix of mapping parameters "F" (from Liu and Kiureghian, 1986) used to generate MVN correlation matrix in \code{CorTransform}.
#' @export
#' @import stats
#' 
#' 
L86 <- function(delta_j, delta_i, rho_ij, distributions){
  
  ## Notation (Liu and Kiureghian, 1986)
  # (row dist) delta_j = sigma_j/mu_j
  # (col dist) delta_j = sigma_i/mu_i
  # rho_ij = i,j correlation
  
  distnames <- c("norm", "unif", "gamma", "lnorm")
  Fmat <- matrix(NA, nrow=4, ncol=4, dimnames = list(distnames,distnames))
  
  # Xj = N and Xi from Group 1
  Fmat["norm","norm"] <- 1
  Fmat["unif","norm"] <- Fmat["norm","unif"] <-1.023
  
  # Xj from Group 2 and Xi = norm
  Fmat["gamma", "norm"] <- 1.001-0.007*delta_j + 0.118*(delta_j^2)
  Fmat["lnorm", "norm"] <- delta_j / sqrt(log(1+(delta_j^2)))
  #
  Fmat["norm", "gamma"] <- 1.001-0.007*delta_i + 0.118*(delta_i^2)
  Fmat["norm", "lnorm"] <- delta_i / sqrt(log(1+(delta_i^2)))
  
  # Xj and Xi from Group 1
  Fmat["unif","unif"] <- 1.047 - 0.047*(rho_ij^2)
  
  # Xj from Group 2, Xi from Group 1
  Fmat["lnorm","unif"] <- 1.019 + 0.0148*delta_j + 0.010*(rho_ij^2) + 0.249*(delta_j^2)
  Fmat["gamma", "unif"] <- 1.023 - 0.007*delta_j + 0.002*rho_ij + 0.127*(delta_j^2)
  #
  Fmat["unif","lnorm"] <- 1.019 + 0.0148*delta_i + 0.010*(rho_ij^2) + 0.249*(delta_i^2)
  Fmat["unif", "gamma"] <- 1.023 - 0.007*delta_i + 0.002*rho_ij + 0.127*(delta_i^2)
  
  # Xj and Xi from Group 2
  if (sum(distributions == "lnorm") >= 2){
    Fmat["lnorm", "lnorm"] <- log(1+(rho_ij*delta_i*delta_j)) / (rho_ij*sqrt(log(1+(delta_i^2))*log(1+(delta_j^2))))
  }
  Fmat["gamma","gamma"] <- 1.002 + 0.022*rho_ij - 0.012*(delta_i + delta_j) + 0.001*(rho_ij^2) + 0.125*((delta_i^2) + (delta_j^2)) - 0.077*rho_ij*(delta_i + delta_j) + 0.014*delta_i*delta_j
  Fmat["gamma", "lnorm"] <- 1.001 + 0.033*rho_ij + 0.004*delta_i - 0.016*delta_j + 0.002*(rho_ij^2) + 0.223*(delta_i^2) + 0.130*(delta_j^2) - 0.104*rho_ij*delta_i + 0.029*delta_i*delta_j - 0.119*rho_ij*delta_i
  #
  Fmat["lnorm", "gamma"] <- 1.001 + 0.033*rho_ij + 0.004*delta_j - 0.016*delta_i + 0.002*(rho_ij^2) + 0.223*(delta_j^2) + 0.130*(delta_i^2) - 0.104*rho_ij*delta_j + 0.029*delta_j*delta_i - 0.119*rho_ij*delta_j
  
  return(Fmat)
}