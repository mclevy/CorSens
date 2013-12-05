#' @title Monte Carlo Estimation of Sensitivity Indices for Models with Dependent Variables (Gaussian Copula Scheme by Kucherenko et al. 2012)
#'
#' @description The \code{K12} function calculates variance decomposition based first order (individual) and interaction (total) sensitivity indices (SI) for models with dependent input variables (data or parameters) using Monte Carlo estimation methods, and a Gaussian copula based approach for sampling from arbitrary multivariate distributions. 
#' 
#' @details The \code{K12} method is a generalization of the method used to compute Sobol' SI, wherein SI respond to the relative impact of interactions and input variable correlations on model output, instead of assuming variable dependence (as in previous Sobol' index estimation methods). \code{K12} relies on specification of a parametric marginal distribution form to each model input variable. Marginal distributions may be any combination of the following distributions: uniform, normal, lognormal, or gamma. Simple random sampling is used for Monte Carlo estimation of SI, and a sample size of no less than \emph{N=10^4} (default) is recommended based on analytical tests. Any number of input variables may be used, however using the recommended sample size with \emph{k>10} input variables may be slow; the function estimates \emph{2k} indices (for first order and total SI) at a cost of \emph{N*(2k+1)} model evaluations. \emph{(Future iterations will implement sampling methods that require fewer samples to achieve SI estimate convergence according to analytical tests.)}
#' 
#' @param data \code{NULL} (default) A matrix with columns taken as model input variables that follow arbitrary marginal distributions (specified by \code{distributions}), and rows as data observations. \code{data} is used only to extract variable means, and variable correlation for transformation to MVN simulated data usign copula based methods. \code{data} is not entered when providing simulation parameters \code{k,sigma,mu} for variables distributions instead.
#' @param distributions A list of distribution names from the list c("unif", "norm", "lnorm", "gamma"), the length of \code{ncol(data)} or \code{k} (in the case of a simulation). Default is all variables are assigned a uniform marginal distribution.
#' @param k \code{NULL} (default) The number of marginal distributions, specified when no data is given, and when \code{sigma} and \code{mu} are provided (such as in a simulation).
#' @param sigma \code{NULL} (default) The marginal distributions' correlation matrix, specified when no data is given, and when \code{k} and \code{mu} are provided (such as in a simulation).
#' @param mu \code{NULL} (default) The marginal distributions' mean vector (of length \code{k}), specified when no data is given, and when \code{k} and \code{sigma} are provided (such as in a simulation).
#' @param model A function passed into \code{K12}. The model must take input variables by name/column number matching the name/column order in \code{data}, or in the order relative to parameters \code{distributions, parameters, sigma, mu} in the case of a simulation.
#' @param parameters \code{NULL} (default) A list containing \code{k} vectors of length two for marginal distribution fit parameters (e.g. \code{mean,sd} for the normal distribution). Parameters are calculated by defualt wtihin \code{K12} using \code{ParamFit} (maximum likelihood fit using \code{fitdistr} from the \code{MASS} package) when \code{data} is supplied. In the case of a simulation, a \code{parameters} list must be supplied.
#' @param N Sample size, by default \emph{=10^4}.
#' @param threshold A vector of any length of desired (arbitrary) sensitivity index thresholds (for return of a \code{$results} matrix indicating variables whose first order and total SI meet the threshold value. By default, the threshold is set to 0.01,0.05, and 0.10 (1\%, 5\%, and 10\% of total output variance). 
#' 
#' @references Kucherenko, S., S. Tarantola, and P. Annoni. "Estimation of Global Sensitivity Indices for Models with Dependent Variables." Computer Physics Communications 183, no. 4 (April 2012): 937-946. doi:10.1016/j.cpc.2011.12.020.
#' @references Saltelli, Andrea, Paola Annoni, Ivano Azzini, Francesca Campolongo, Marco Ratto, and Stefano Tarantola. "Variance Based Sensitivity Analysis of Model Output. Design and Estimator for the Total Sensitivity Index." Computer Physics Communications 181, no. 2 (February 2010): 259-270. doi:10.1016/j.cpc.2009.09.018.
#' @references Saltelli, Andrea. "Making Best Use of Model Evaluations to Compute Sensitivity Indices." Computer Physics Communications 145, no. 2 (May 15, 2002): 280-297. doi:10.1016/S0010-4655(02)00280-1.
#' @references Sobol', I.M. "Global Sensitivity Indices for Nonlinear Mathematical Models and Their Monte Carlo Estimates." Mathematics and Computers in Simulation 55, no. 1-3 (February 15, 2001): 271-280. doi:10.1016/S0378-4754(00)00270-6.
#' 
#' @return list \code{K12} returns a list containing the first order SI estimates \code{$SY}; the total SI estimates \code{$STy}; and a printed matrix \code{$results} of first order and total SI achieving threshold values.
#' @export
#' @importFrom MASS fitdistr
#' @import stats
#' @import matrixcalc
#' 
#' @examples \dontrun{
#' # Data (trees) Example
#' cor(trees) # dependent variables
#' mod_trees <- function(data) data[ ,1]*data[ ,2]*data[ ,3]
#' s_trees <- K12(data = trees, distributions = c("gamma", "gamma", "norm"), model = mod_trees)
#' s_trees$SY
#' s_trees$STy
#' s_trees$results
#' 
#' # Simulation Example
#' N <- 10^4 # observations
#' k <- 5 # variables
#' sigma <- diag(k) # independent
#' sim <- matrix(runif(N*k, min=0, max=10), ncol=k)
#' mod_sim <- function(data) (data[ ,1]*data[ ,2]) + data[ ,3] + data[ ,4] + data[ ,5]
#' s_sim <- K12(data = sim, model = mod_sim) # distributions = uniforms by default
#' s_sim$SY
#' s_sim$STy
#' s_sim$results
#' }
#' 
#' 
K12 <- function(data=NULL, distributions=NULL, parameters = NULL, N = 10^4, sigma = NULL, model, k=NULL, mu = NULL, threshold = c(0.01,0.05,0.1)){
  
  if(is.null(data) && is.null(sigma)) stop("If no data is provided, must specify sigma matrix")
  if(is.null(data) == T && is.null(parameters)) stop("If no data is provided, must specify parameters for chosen (default uniform) distributions")
  if(is.null(data) == T && is.null(mu)) stop("If no data is provided, must specify mean (vector) for distributions")
  if(is.null(data) == T && is.null(k)) stop("If no data is provided, must specify k = the number of data variables")
  if(is.null(model)) stop("Model function is missing")
  
  if (is.null(sigma) == F){
    # convert to correlation; if already cor, doesn't matter
    sigma <- cov2cor(sigma)
  } else if (is.null(sigma) == T){
    sigma <- cor(data,use="pairwise.complete.obs")
  }
  
  if(is.positive.semi.definite(round(sigma,2)) == F) stop("Correlation matrix is not positive semi-definite.")
  
  if (is.null(k) == T){
    k <- dim(data)[2]
    test <- F
  } else {
    k <- k
    test <- T
  }
  
  q <- k-1
  
  ## specify distributions according to inputs (data or mu,sigma,k specified for a simulated distribution)
  
  if (is.null(distributions) && (is.null(data) == F)){
    distributions <- rep("unif", ncol(data))
  } else if (is.null(distributions) && (is.null(data) == T)){
    distributions <- rep("unif", k)
  } else {
    distributions <- distributions
  }
  
  if (is.null(data) == T){  # no data matrix
    if(k != length(distributions)) stop("Number of variables (k) does not match the length of the list of distributions") 
  } else if (is.null(data) == F){
    if(ncol(data) != length(distributions)) stop("Number of variables (data columns) does not match the length of the list of distributions") 
  }
  
  ##  extract arbitrary distribution parameters for use in later transformations
  
  if (is.null(parameters) == T){
    params <- ParamFit(data, distributions)
  } else {
    params <- parameters
  }
  
  ## Specify transformed sigma; from data and distributions, or entered manually
  if (is.null(sigma) == T){
    sigma <- CorTransform(data, k=k, distributions)
  } else if (is.null(data) == F && is.null(sigma) == F){
    sigma <- CorTransform(data=data, k=k, distributions, sigma)
  } else if (is.null(data) == T && is.null(sigma) == F && is.null(mu) == F){
    sigma <- CorTransform(distributions=distributions,k=k, sigma=sigma, mu=mu)
  }
  
  #########################################################################
  ### Set up for: Construct sets (y,z), (y, z_bp) and complimentary set (y_bp, z)
  ## x = (y,z), x_p = (y_tp, z_tp) ; t=tilde, p=prime (paper notation)
  ## x= (y, z) with rotating col=i as y; the rest (-i) are z.
  
  # x matrix contains CORRELATED variables; gen normals and transform using sigma
  u <- matrix(runif((k) * N), nrow = N) # k-2-dim (k-2 = correlated vars only) ; uniform sample, u
  xt <- qnorm(u) # normal inverse transform
  
  if (is.null(mu) == T){
    mu_x <- c(colMeans(data, na.rm=T))
  }else{
    mu_x <- mu
  }
  
  A <- t(chol(sigma)) # check: A %*% t(A) = sigma ; sigma is just correlated vars
  x <- t( apply(A %*% t(xt), MARGIN=2, FUN=function(x){x + mu_x}))
  # pairs(x[sample(nrow(x),10^3), ])
  
  ## x_p = cbind(y_tp, z_tp) with rotating col=i as y_tp; the rest (-i) cols are z_tp
  
  # x_tp matrix is UNCORRELATED standard normal variables = no transformation using sigma
  u_p <- matrix(runif(k * N), nrow = N) # k-dim uniform sample, u_p
  xt_p <- qnorm(u_p) # transform to standard normal
  # pairs(xt_p[sample(nrow(x),10^3), ])
  
  ## Calculate mu's = mu_y, mu_z: For matrix x, mu_y is colMeans of the y "subset" col (i) for each subset; mu_z is the colMeans of the remaining (-i) cols of x.
  
  Mu_x <- colMeans(x)
  # Mu_xt_p <- colMeans(xt_p) # not used below
  
  #########################################################################
  ### Construct set (y, z_bp)
  
  ## compute mu_zc using y; from eqn.(3.4): mu_zc = mu_z + sigma_zy %*% solve(sigma_y) %*% (y-mu_y)
  
  st <- seq(1,q*k,by=q) ; end <- seq(q,k*q, by=q) # index
  Mu_zc <- matrix(NA, nrow = q*k, ncol=N)
  # computes Mu_zc for each subset i = 1:3, for each N samples
  for (i in 1:k){
    # partitioned sigma matrix: same for all subsets of x
    sigma_y <- matrix(sigma[i,i])
    sigma_zy <- matrix(sigma[i,-i])
    Mu_zc[st[i]:end[i], ] <- apply(sigma_zy %*% solve(sigma_y) %*% (x[ ,i] - Mu_x[i]), MARGIN=2, FUN=function(x){x + matrix(Mu_x[-i])})
  }
  
  ##  compute z_bp (bp = bar, prime); zt_p = xt_p[ ,-i]; z_bp = A_zc %*% z_tp + mu_zc; which follows (q = k-s)-dim conditional normal distribution in eq.(3.3)
  
  Z_bp <- matrix(NA, nrow = q*k, ncol=N)
  for (i in 1:k){
    sigma_y <- sigma[i,i]
    sigma_yz <- sigma[-i,i]
    sigma_zy <- sigma[i,-i]
    sigma_z <- sigma[-i,-i]
    sigma_zc <- sigma_z - sigma_zy %*% solve(sigma_y) %*% sigma_yz
    A_zc <- t(chol(sigma_zc))
    Z_bp[st[i]:end[i], ] <- apply(xt_p[ ,-i], MARGIN=c(1), FUN=function(x){A_zc %*% x}) + Mu_zc[st[i]:end[i], ]
  }
  
  ## create vector (y,z_bp)
  
  nst <- seq(1,k*N,by=N) ; nend <- seq(N,k*N, by=N) 
  
  Y.Z_bp <- matrix(NA, nrow = k*N, ncol=k)
  for (i in 1:k){
    Y.Z_bp[nst[i]:nend[i],i] <- x[,i]
    Y.Z_bp[nst[i]:nend[i],-i] <- t(Z_bp[st[i]:end[i], ])
  }
  
  ## copula: transform normals to arbitrary marginal distribution
  # Normal CDF
  Y.Z_bp <- Fnorm(Y.Z_bp) # Uniforms
  
  # Arbitrary quantile
  Y.Z_bp <- Ginv(Y.Z_bp, distributions, params)
  
  # Y.Z_bp: subsets of 1:N each, with rotating col i as the original col from x, z from Z_bp transform of z (-i) cols of xt_p
  
  #########################################################################
  ### Construct (complimentary) set (y_bp, z)
  ## From above: xt_p = cbind(y_tp, z_tp) with rotating col=i as y_tp => y_tp = xt_p[ ,i]; the rest (-i) are z_tp
  ## compute mu_yc using z (from x[ ,-i]); mu_yc = mu_y + sigma_yz %*% solve(sigma_z) %*% (z-mu_z)
  
  Mu_yc <- matrix(NA, nrow = 1*k, ncol=N)
  for (i in 1:k){
    sigma_yz <- sigma[-i,i]
    sigma_z <- sigma[-i,-i]
    Mu_yc[i, ] <- apply(sigma_yz %*% solve(sigma_z) %*% apply(x[ ,-i], MARGIN=1, FUN=function(x){x - Mu_x[-i]}), MARGIN=2, FUN=function(x){x + Mu_x[i]})
  }
  
  ## compute y_bp (bp = bar, prime); y_bp = A_yc %*% y_tp + mu_yc; which follows 1-dim conditional normal distribution
  
  Y_bp <- matrix(NA, nrow = 1*k, ncol=N)
  for (i in 1:k){
    # partitioned sigma matrix: same for all subsets of x
    sigma_y <- sigma[i,i]
    sigma_yz <- sigma[-i,i]
    sigma_zy <- sigma[i,-i]
    sigma_z <- sigma[-i,-i]
    sigma_yc <- sigma_y - sigma_yz %*% solve(sigma_z) %*% sigma_zy
    A_yc <- t(chol(sigma_yc)) # A_yc %*% t(A_yc)
    Y_bp[i, ] <- apply(matrix(xt_p[ ,i]), MARGIN=c(1), FUN=function(x){A_yc %*% x} ) + Mu_yc[i, ]
  }
  
  ##  create vector cbind(y_bp,z)
  
  Y_bp.Z <- matrix(NA, nrow = k*N, ncol=k)
  for (i in 1:k){
    Y_bp.Z[nst[i]:nend[i],i] <- Y_bp[i, ]
    Y_bp.Z[nst[i]:nend[i],-i] <-  x[,-i]
  }
  
  ## copula: transform normals to arbitrary marginal distribution
  # Normal CDF
  Y_bp.Z <- Fnorm(Y_bp.Z) # Uniforms

  # Arbitrary quantile
  Y_bp.Z <- Ginv(Y_bp.Z, distributions, params)
  
  # Y_bp.Z: subsets of 1:N each, with rotating col i as the Y_bp col from transformed xt_p (i) cols; z from original x (-i)
  
  #########################################################################
  ### Construct joint set (y, z)
  # Joint set (y,z) is just the original x matrix x=(y,z); the "y" will be the rotating col i, and the z the (-i) cols for i = 1,...,k
  
  ## copula: transform normals to arbitrary marginal distribution
  # Normal CDF
  Y.Z <- Fnorm(x) # Uniforms
  # Arbitrary quantile
  Y.Z <- Ginv(Y.Z, distributions, params)
  
  #########################################################################
  ### Calculate Sensitivity Indices SY and STy
  
  mod <- model
  
  f.yz <- matrix(mod(Y.Z)) # ouput matrix (N x 1)
  f.y.z_bp <- matrix(mod(Y.Z_bp), ncol=k) # output matrix (N x k) for each k subset
  f.y_bp.z <- matrix(mod(Y_bp.Z), ncol=k) # output matrix (N x k) for each k subset
  
  f0 <- sum(f.yz)/N
  D <- sum(f.yz^2)/N - f0^2
  
  Sy <- apply(f.y.z_bp, MARGIN=2, FUN=function(x){(sum(f.yz * x)/N - f0^2)/D})
  STy <- apply(f.y_bp.z, MARGIN=2, FUN=function(x){(sum((f.yz - x)^2)/N)/(2*D)})
  
  # results table with thresholds
  kk <- list(Sy, STy)
  names(kk)=c("SY", "STy")
  
  thresh <- threshold
  threshnames_i <- paste(thresh, "_i", sep="")
  threshnames_T <- paste(thresh, "_T", sep="")
  results_names <- c("Rank_i", threshnames_i, "Rank_T", threshnames_T)
  vi <- seq(1:k)
  oi <- order(kk$SY, decreasing = T)
  ot <- order(kk$STy, decreasing = T)
  
  if (test == T){
    nms <- as.character(1:k)
  } else {
    nms <- names(data)
  }
  
  results <- matrix(0, nrow = k, ncol=length(results_names), dimnames = list(nms, results_names))
  results <- as.data.frame(results)
  
  for (i in 1:k){
    results[oi[i],"Rank_i"] <- vi[i]
    results[ot[i],"Rank_T"] <- vi[i]
  }
  
  for (i in 1:length(thresh)){
    if (sum(kk$SY > threshnames_i[i]) > 0){
      results[kk$SY > threshnames_i[i], threshnames_i[i]] <- "X"
    }
    if (sum(kk$STy > threshnames_T[i]) > 0){
      results[kk$STy > threshnames_T[i], threshnames_T[i]] <- "X"
    }
  }
  
  results[results == 0] <- "."
  
  out=list(Sy, STy, results)
  names(out)=c("SY", "STy", "results")
  
  return(out)
  
}