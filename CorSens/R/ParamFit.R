#' @title Extract distribution parameter fits to empirical PDFs
#'
#' @description The function extracts the Maximum Likelihood (ML) parameter fits to empirical, marginal data distributions from the normal, lognormal, and gamma distributions using the \code{fitdistr} function from the \code{MASS} package, and the simple minimum and maximum for data assigned a uniform distribution.
#' 
#' @details The function is used internal to the \code{K12} function or may be used independently.
#'
#' @param data A data matrix with columns taken as arbitrary marginal distributions to be fit with parameters specific to distributions specified by the \code{distributions} vector.
#' @param distributions A vector of distribution names from the list c("unif", "norm", "lnorm", "gamma"), the length of ncol(data).
#'
#' @return list A list of arbitrary distribution fit parameters
#' @export
#' @importFrom MASS fitdistr
#' @import stats
#' 
#' 
ParamFit <- function (data, distributions){
  
  if(sum(is.na(charmatch(distributions, c("unif", "norm", "lnorm", "gamma")))) != 0) stop("Distribution names do not match available distributions") 
  
  dist <- distributions
  
  # re-name distributons to match fitdistr (MASS) distribution names
  dist[dist == "norm"] <- "normal"
  dist[dist == "lnorm"] <- "lognormal"
  dist[dist == "gamma"] <- "gamma"
  params <- vector("list", length(dist))
  
  for (i in 1:length(distributions)){
    if (dist[i] != "unif" && dist[i] != "gamma"){
      params[[i]] <- as.numeric(c(fitdistr(data[ ,i], dist[i]))$estimate)
    } else if (dist[i] == "gamma"){
      params[[i]] <- as.numeric(c(fitdistr(data[ ,i], dist[i], lower = 0.001))$estimate)
    } else if (dist[i] == "unif"){
      params[[i]] <- c(min(data[ ,i], na.rm=T), max(data[ ,i], na.rm=T)) # manually fit uniform with naive param est
    }
    
  }
  return(params) # a list
}