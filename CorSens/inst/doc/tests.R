
## ----, echo = FALSE------------------------------------------------------

rm(list= ls());
library(CorSens)

N <- 10^5


## ----, echo = FALSE------------------------------------------------------

############################ Test Case 1

# dimensions
k <- 3 # dimensions (variables/parameters) # test case 1
q <- k-1 # remaining vars from each subset

# model
mod1 <- function(data){
  Y <- (data[ ,1] + data[ ,2] + data[ ,3])
  return(Y)
}

## sigma matrices for different rho values
sig1 <- function(rho, mu, std, k=3){ # cor, mean, sd, # of vars
  sigma <- matrix(rep(0, k^2), nrow=k) # cor mat init
  diag(sigma) <- 1
  sigma[k,k] <- std^2
  sigma[3,2] <- rho* std ; sigma[2,3] <- rho*std
  return(sigma)
}

rhos1 <- c(0,0.5,-0.5,0.8,-0.8)
mus <- rep(0,k)
std <- 2

Ki <- KT <- KAi <- KAT <-  NA

for (i in 1:length(rhos1)){
  
  # Sigma matrix
  sigma <- sig1(rhos1[i],0,std)
  
  p <- list(c(mus[1],sqrt(sigma[1,1])), c(mus[2],sqrt(sigma[2,2])), 
            c(mus[3],sqrt(sigma[3,3])))
  
  K <- K12(N = N, distributions = rep("norm",k), sigma = sigma, 
           model=mod1, k = 3, mu = mus, parameters=p)
  
  Ki <- c(Ki,K$SY)
  KT <- c(KT,K$STy)
  
  ## Analytical values for linear model f(x1 + x2 + x3) with sigma above
  s1 <- 1 / (2 + std^2 + 2*rhos1[i]*std)
  s2 <- (1 + rhos1[i]*std)^2 / (2 + std^2 + 2*rhos1[i]*std)
  s3 <- (std + rhos1[i])^2 / (2 + std^2 + 2*rhos1[i]*std)
  KAi <- c(KAi,s1,s2,s3)
  
  s1T <- 1 / (2 + std^2 + 2*rhos1[i]*std)
  s2T <- (1 - rhos1[i]^2) / (2 + std^2 + 2*rhos1[i]*std)
  s3T <- (std^2) * (1 - rhos1[i]^2) / (2 + std^2 + 2*rhos1[i]*std)
  KAT <- c(KAT,s1T,s2T,s3T)
  
  # print(i)
}

t1 <- data.frame(variable = rep(c("x1", "x2", "x3"),5), 
                 rho = rep(rhos1, each=3, 1), Si=NA, Si.A=NA, 
                 STi=NA, STi.A=NA)

# enter Kucherenko MC estimates
t1[ ,3] <- Ki[-1]
t1[ ,5] <- KT[-1]
t1[ ,4] <- KAi[-1]
t1[ ,6] <- KAT[-1]

# Print table == Table 1 (p.943) in Kucherenko et al. 2012
numVars <- sapply(t1, is.numeric) 
t1[numVars] <- lapply(t1[numVars], round, digits = 2) 
print(t1)



## ----, echo=FALSE--------------------------------------------------------

k <- 4 # dimensions (variables/parameters) # test case 2

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

p <- list(c(mus[1],sqrt(sigma[1,1])), c(mus[2],sqrt(sigma[2,2])), c(mus[3],sqrt(sigma[3,3])), c(mus[4],sqrt(sigma[4,4])) )

t2 <- data.frame(variable = c("x1", "x2", "x3", "x4"), Si=NA, Si.A=NA, STi=NA, STi.A=NA)

K <- K12(N = N, distributions = rep("norm",k), sigma = sigma, model=mod2, k = k, mu = mus, parameters=p)

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

# Print table == Table 2 (p.945) in Kucherenko et al. 2012
numVars <- sapply(t2, is.numeric) 
t2[numVars] <- lapply(t2[numVars], round, digits = 2) 
print(t2)



## ----, echo=FALSE, fig.width=7, fig.height=5-----------------------------

############################ Test Case 3

k <- 3

# data
u <- matrix(runif(N*k, min=-pi, max = pi), nrow=N)
mu_x <- colMeans(u)

# correlation
corseq <- seq(-0.9, 0.9,length.out=10)

#toy model
ishigami <- function(data){
  Y <- ( sin(data[ ,1]) + (7*(sin(data[ ,2]))^2) 
         + (0.1* ((data[ ,3])^4) * sin(data[ ,1])) )
  return(Y)
}

# initialize
ki_mat <- matrix(NA, nrow = k, ncol=length(corseq), 
                 dimnames = list(1:k, corseq)) # indiv index matrix
kt_mat <- matrix(NA, nrow = k, ncol=length(corseq), 
                 dimnames = list(1:k, corseq)) # total index matrix

for (i in 1:length(corseq)){
  
  # Desired correlation (of input data variables) = actual cor
  rho <- 0.0
  sigma <- matrix(rep(NA, k^2), nrow=k) # cor mat init
  std <- 1
  diag(sigma) <- std^2 ; sigma[lower.tri(sigma)] <- rho ; 
  sigma[upper.tri(sigma)] <- rho
  sigma[1,3] <- sigma[3,1] <- corseq[i]
  
  ki <- suppressWarnings(K12(data = u, N=N, model=ishigami, sigma = sigma))
  # suppress warnings because sigma has all(diag(sigma) == 1), which triggers warning message about checking that cov matrix supplied (instead of cor)
  
  ki_mat[ ,i] <- ki$SY
  kt_mat[ ,i] <- ki$STy
  # print(i)
  
}


## plot

plot(corseq, ki_mat[1, ], ylab="Sensitivity Index Value", type="b", 
     pch=21, lty=1, ylim=c(0,0.7), col=1, xlab=expression(rho["1,3"]))

lines(corseq, kt_mat[1, ], ylab="Index Value", type="b", 
      pch=21, lty=2, col=1)
lines(corseq, ki_mat[2, ], ylab="Index Value", type="b", 
      pch=24, lty=1, col=2)
lines(corseq, kt_mat[2, ], ylab="Index Value", type="b", 
      pch=24, lty=2, col=2)
lines(corseq, ki_mat[3, ], ylab="Index Value", type="b", 
      pch=23, lty=1, col=4)
lines(corseq, kt_mat[3, ], ylab="Index Value", type="b", 
      pch=23, lty=2, col=4)

legend("top", expression("x1"["i"],"x1"["T"], "x2"["i"],"x2"["T"], 
       "x3"["i"],"x3"["T"]), pch=c(21,21,24,24,23,23), 
       lty = c(1,2,1,2,1,2), col=c(1,1,2,2,4,4), pt.bg = "white", 
       bg="white", horiz=T, bty="n", ncol=1, y.intersp = 1)



