context("This performs tests on CorSens functions")

test_that("Ensuring that CorSens functions: K12 and functions operating internal to K12, produce desired results based on analytical tests from Kucherenko et al. 2012",{
  
  # Set Up: Test Case 2 from Kucherenko et al. 2012
  N <- 10^4
  k <- 4
  
  # model
  mod2 <- function(data){
    Y <- (data[ ,1]*data[ ,3] + data[ ,2]*data[ ,4])
    return(Y)
  }
  
  # Sigma matrix
  sigma <- matrix(rep(0, k^2), nrow=k) # cor mat init
  diag(sigma) <- c(16,4,4*10^4, 9*10^4)
  sigma[1,2] <- sigma[2,1] <- 2.4 ;
  sigma[3,4] <- sigma[4,3] <- -1.8*(10^4);
  mus <- c(0,0,250,400)
  
  p <- list(c(mus[1],sqrt(sigma[1,1])), c(mus[2],sqrt(sigma[2,2])), 
            c(mus[3],sqrt(sigma[3,3])), c(mus[4],sqrt(sigma[4,4])) )
  
  t2 <- data.frame(variable = c("x1", "x2", "x3", "x4"), Si=NA, 
                   Si.A=NA, STi=NA, STi.A=NA)
  
  K <- K12(N = N, distributions = rep("norm",k), sigma = sigma, 
           model=mod2, k = k, mu = mus, parameters=p)
  
  # Analytical results:
  rho12 <- sigma[1,2]/sqrt(sigma[1,1]*sigma[2,2])
  rho34 <- sigma[3,4]/sqrt(sigma[3,3]*sigma[4,4])
  D <- sigma[1,1]*(sigma[3,3] + mus[3]^2) + sigma[2,2]*(sigma[4,4] + mus[4]^2) + 2*sigma[1,2]*(sigma[3,4] + (mus[3]*mus[4]))
  
  s1 <- (sigma[1,1]*(mus[3] + mus[4]*rho12*sqrt(sigma[2,2]/sigma[1,1]))^2)/D
  s2 <- (sigma[2,2]*(mus[4] + mus[3]*rho12*sqrt(sigma[1,1]/sigma[2,2]))^2)/D
  s3 <- 0
  s4 <- 0
  
  s1T <- (sigma[1,1]*(1-rho12^2)*(sigma[3,3] + mus[3]^2))/D
  s2T <- (sigma[2,2]*(1-rho12^2)*(sigma[4,4] + mus[4]^2))/D
  s3T <- (sigma[1,1]*sigma[3,3]*(1-rho34^2))/D
  s4T <- (sigma[2,2]*sigma[4,4]*(1-rho34^2))/D
  
  # enter Kucherenko MC and Analytical estimates
  t2[ ,2] <- K$SY
  t2[ ,3] <- c(s1,s2,s3,s4)
  t2[ ,4] <- K$STy
  t2[ ,5] <- c(s1T,s2T,s3T,s4T)
  
  # Table 2 (p.945) in Kucherenko et al. 2012
  numVars <- sapply(t2, is.numeric) 
  t2[numVars] <- lapply(t2[numVars], round, digits = 1) 
  
  # Test match of estimated and analytical values
  expect_that(t2[,2], equals(t2[,3]))
  expect_that(t2[,4], equals(t2[,5]))
  
})
