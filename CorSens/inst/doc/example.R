## ----, echo=TRUE---------------------------------------------------------

set.seed(1208)

# dimensions of data
N <- 10^5
k <- 4

## Simulate correlated data:

# Generate normals
u <- matrix(runif((k) * N), nrow = N) # uniforms
xt <- qnorm(u) # normal inverse transform
mu_x <- rep(0,4) # mean vector

# Correlation matrix (zero correlation)
rho <- 0.0
sigma <- matrix(rep(NA, k^2), nrow=k)
std <- 1
diag(sigma) <- std^2 ; sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)] <- rho

# Generate uniform, correlated variables with cholesky decomposition
A <- t(chol(sigma))
x <- t( apply(A %*% t(xt), MARGIN=2, FUN=function(x){x + mu_x}))
U <- apply(x,2,pnorm)

# simulate data according to selected arbitrary distributions
x1 <-  qgamma(U[,1], shape=5, rate=1/2) # gamma = exponential
x2 <- qnorm(U[,2]) # normal
x3 <- rlnorm(N, 0, mean(qlnorm(U[,3]))*0.2) # lognormal (meets allowable coefficient of correlation range sigma/mu = [0.1,0.5] accroding to Liu et. al., 1986)
x4 <- qunif(U[,4]) # uniform

# matrix of simulated, correlated data (4 variables)
data <- cbind(x1,x2,x3,x4)

## ----, echo=FALSE, fig.width=7, fig.height=5-----------------------------
par(mfrow=c(2,2))

for (i in 1:ncol(data)){
  hist(data[ ,i], freq=F, main=NULL, xlab=colnames(data)[i])
  lines(density(data[ ,i]), col=2, lwd = 2)
}

## ----, echo=FALSE--------------------------------------------------------
round(cov(data),1)

## ----, echo=TRUE---------------------------------------------------------
# Assign distribution manually
distributions <- c("gamma", "norm","lnorm", "unif")

# Toy model
toy <- function(data){
  Y <- (data[ ,1]*data[ ,2]) + (data[ ,1]*data[ ,4])^2 + data[ ,1] + data[ ,2] + data[ ,3] + data[ ,4]
  return(Y)
}

## ----, echo=TRUE---------------------------------------------------------
### Sobol' 2002 method
library(sensitivity, quietly=T)

# this function requires making two separate samples (a,b)
n_samps <- N
a <- data[sample(nrow(data),size=n_samps,replace=TRUE), ]
b <- data[sample(nrow(data),size=n_samps,replace=TRUE), ]
s1_ex <- sobol2002(model = toy, X1 = a, X2 = b)

# Individual indices
s1_ex$S
# Total indices
s1_ex$T

## ----, echo=TRUE---------------------------------------------------------
### Kucherenko 2012 method
library(CorSens)
k1_ex <- K12(N=N, data = data, distributions=distributions, model = toy)

# Individual indices
k1_ex$SY
# Total indices
k1_ex$STy
# sensitivity according to (arbitrary) thresholds of significance (defaults)
k1_ex$results

## ----, echo=FALSE--------------------------------------------------------
### Plot distributions function

ex_plot <- function(k,s){
  matplot(cbind(k$SY,s$S[1]), xaxt = "n", ylab="Index Value", pch=c(21,22), main="Individual Indices", ylim=c(0,1))
  axis(1, 1:4, c("x1", "x2", "x3", "x4"))
  abline(v = 1:4, col = "lightgray", lty=3)
  legend("topright", c("Kucherenko", "Sobol'"), pch=c(21,22), col=1:2, cex=1.2)
  
  matplot(cbind(k$STy,s$T[1]), xaxt = "n", ylab="Index Value", pch=c(21,22), main="Total Indices", ylim=c(0,1))
  axis(1, 1:4, c("x1", "x2", "x3", "x4"))
  abline(v = 1:4, col = "lightgray", lty=3)
  legend("topright", c("Kucherenko", "Sobol'"), pch=c(21,22), col=1:2, cex=1.2)
}

## ----, echo=FALSE, fig.width=7, fig.height=4-----------------------------
par(mfrow=c(1,2))
ex_plot(k1_ex, s1_ex)

## ----, echo=TRUE---------------------------------------------------------
# Correlation matrix (WITH correlation)
rho <- 0.0
sigma <- matrix(rep(NA, k^2), nrow=k)
std <- 1
diag(sigma) <- std^2 ; sigma[lower.tri(sigma)] <- rho ; 
sigma[upper.tri(sigma)] <- rho

rhos1 <- 0.8
sigma[2,1] <- sigma[1,2] <- rhos1

A <- t(chol(sigma))
x <- t( apply(A %*% t(xt), MARGIN=2, FUN=function(x){x + mu_x}))
U <- apply(x,2,pnorm) # Uniform, correlated variables

x1 <-  qgamma(U[,1], shape=5, rate=1/2) # gamma
x2 <- qnorm(U[,2]) # normal
x3 <-  rlnorm(N, 0, mean(qlnorm(U[,3]))*0.2) # lognormal (meets allowable coefficient of correlation range sigma/mu = [0.1,0.5] accroding to Liu et. al., 1986)
x4 <- qunif(U[,4]) # uniform

data <- cbind(x1,x2,x3,x4) # simulated, correlated data

round(cov(data),1) # covariance matrix

## ----, echo=FALSE, fig.width=7 ,fig.height=4-----------------------------
### Kucherenko 2012 method
k2_ex <- K12(N=N, data = data, distributions=distributions, model = toy)

### Sobol' 2002 method
a <- data[sample(nrow(data),size=n_samps,replace=TRUE), ]
b <- data[sample(nrow(data),size=n_samps,replace=TRUE), ]
s2_ex <- sobol2002(model = toy, X1 = a, X2 = b)

### Plot Comparison

par(mfrow=c(1,2))
ex_plot(k2_ex, s2_ex)


## ----, echo=TRUE---------------------------------------------------------
k2_ex$results

## ----, echo=TRUE---------------------------------------------------------
############ Range of Correlation Example

# correlation
corseq <- seq(-0.9, 0.9,length.out=10)

# initialize
ki_mat <- matrix(NA, nrow = k, ncol=length(corseq), 
                 dimnames = list(1:k, corseq)) # indiv index matrix
kt_mat <- matrix(NA, nrow = k, ncol=length(corseq), 
                 dimnames = list(1:k, corseq)) # total index matrix

for (i in 1:length(corseq)){
  
  # Correlation matrix (WITH CORRELATION)
  rho <- 0.0
  sigma <- matrix(rep(NA, k^2), nrow=k) # cor mat init
  std <- 1
  diag(sigma) <- std^2 ; sigma[lower.tri(sigma)] <- rho ; 
  sigma[upper.tri(sigma)] <- rho
  
  sigma[2,1] <- sigma[1,2] <- corseq[i]
  
  A <- t(chol(sigma))
  x <- t( apply(A %*% t(xt), MARGIN=2, FUN=function(x){x + mu_x}))
  
  U <- apply(x,2,pnorm) # Uniform, correlated variables
  
  x1 <-  qgamma(U[,1], shape=5, rate=1/2) # gamma
  x2 <- qnorm(U[,2]) # normal
  x3 <-  rlnorm(N, 0, mean(qlnorm(U[,3]))*0.2) # lognormal (meets allowable coefficient of correlation range sigma/mu = [0.1,0.5] accroding to Liu et. al., 1986)
  x4 <- qunif(U[,4]) # uniform
  
  data <- cbind(x1,x2,x3,x4) # simulated, correlated real data
  
  ki <- K12(N=N, data = data, distributions=distributions, 
            model = toy)
  
  ki_mat[ ,i] <- ki$SY
  kt_mat[ ,i] <- ki$STy
  # print(i)
  
}


## ----, echo=FALSE, fig.width=7, fig.height=5-----------------------------

par(mfrow=c(1,1))
cols <- c(1,2,3,4) # plot colors
pchs <- c(21,24,23,22) # point types

plot(corseq, ki_mat[1, ], ylab="Sensitivity Index Value", type="b", 
     pch=21, lty=1, ylim=c(0,1), col=1, xlab=expression(rho["x1,x2"]))
lines(corseq, kt_mat[1, ], ylab="Index Value", type="b", 
      pch=21, lty=2, col=1)

for (i in 2:k){
  lines(corseq, ki_mat[i, ], ylab="Index Value", type="b", 
        pch=pchs[i], lty=1, col=cols[i])
  lines(corseq, kt_mat[i, ], ylab="Index Value", type="b", 
        pch=pchs[i], lty=2, col=cols[i])
}

legend("top", expression("x1"["i"],"x1"["T"], "x2"["i"],
       "x2"["T"], "x3"["i"],"x3"["T"], "x4"["i"],"x4"["T"]), 
       pch=rep(pchs,each = 2), lty = rep(c(1,2),4), 
       col=rep(cols,each = 2), pt.bg = "white", bg="white", 
       horiz=T, bty="n", ncol=1, y.intersp = 1)



