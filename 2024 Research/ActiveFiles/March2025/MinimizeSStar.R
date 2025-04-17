


x <- runif(1000, min=1, max=3)
x0 = 1.5
bandwidth = 0.4

gaussKernel <- function(x){
  1/sqrt(2*pi) * exp(-x^2 / 2)
}

Snj <- function(j){
  
  output = sum((x-x0)**j * gaussKernel((x-x0)/bandwidth))
  return(output)
}

Sn0 <- Snj(0)
Sn1 <- Snj(1)
Sn2 <- Snj(2)
Sn3 <- Snj(3)
Sn4 <- Snj(4)
Sn5 <- Snj(5)
Sn6 <- Snj(6)
Sn7 <- Snj(7)

sStar <- matrix(c(Sn0, Sn1, Sn3, Sn1, Sn2, Sn4, Sn3, Sn4, Sn6), nrow=3)
sStarInv <- solve(sStar)

S <- matrix(c(Sn0, Sn1, Sn2, Sn3, Sn1, Sn2, Sn3, Sn4, Sn2, Sn3, Sn4, Sn5, Sn3, Sn4, Sn5, Sn6), nrow=4)
sInv <- solve(S)

sStarVector <- matrix(c(Sn4, Sn5, Sn7), ncol=1)
sVector <- matrix(c(Sn4, Sn5, Sn6, Sn7), ncol=1)

bias <- sInv %*% sVector
biasStar <- sStarInv %*% sStarVector

print(bias)
print(biasStar)

runs <- 1000
biasStorage <- matrix(nrow=runs, ncol=4)
biasStarStorage <- matrix(nrow=runs, ncol=3)
for(i in 1:runs){
  
  x <- runif(100)
  x0 = 0.5
  bandwidth = 0.1
  
  
  Sn0 <- Snj(0); Sn1 <- Snj(1); Sn2 <- Snj(2); Sn3 <- Snj(3); Sn4 <- Snj(4)
  Sn5 <- Snj(5); Sn6 <- Snj(6); Sn7 <- Snj(7);
  
  sStar <- matrix(c(Sn0, Sn1, Sn3, Sn1, Sn2, Sn4, Sn3, Sn4, Sn6), nrow=3)
  sStarInv <- solve(sStar)
  
  S <- matrix(c(Sn0, Sn1, Sn2, Sn3, Sn1, Sn2, Sn3, Sn4, Sn2, Sn3, Sn4, Sn5, Sn3, Sn4, Sn5, Sn6), nrow=4)
  sInv <- solve(S)
  
  sStarVector <- matrix(c(Sn4, Sn5, Sn7), ncol=1)
  sVector <- matrix(c(Sn4, Sn5, Sn6, Sn7), ncol=1)
  
  bias <- sInv %*% sVector
  biasStar <- sStarInv %*% sStarVector
  
  biasStorage[i,] <- bias
  biasStarStorage[i,] <- biasStar
}

#SnInverse for Beta0
beta <- 1
par(mfrow=c(2,1))
starHist <- hist(abs(biasStarStorage[,2]))
Hist <- hist(abs(biasStorage[,2]))
par(mfrow=c(1,1))
plot( starHist, col=rgb(0,0,1,1/4), xlim=c(0.0001,0.0005))  # first histogram
plot( Hist, col=rgb(1,0,0,1/4), xlim=c(0.0001,0.0005), add=T)  # sec
     

ts.plot(abs(biasStarStorage[1:50,1]), lwd=2)
lines(abs(biasStorage[1:50,1]), col='red', lwd=2)

print(mean(abs(biasStarStorage[,1])))
print(mean(abs(biasStorage[,1])))



XStar <- matrix(nrow=100, ncol=3)
XStar[,1] <- 1
XStar[,2] <- x - x0
XStar[,3] <- (x-x0)**3

X <- matrix(nrow=100, ncol=4)
X[,1] <- 1
X[,2] <- x - x0
X[,3] <- (x-x0)**2
X[,4] <- (x-x0)**3

W <- diag(gaussK((x-x0)/bandwidth))


variance <- sInv %*% t(X) %*% W %*% W %*% X %*% sInv
varianceStar <- sStarInv %*% t(XStar) %*% W %*% W %*% XStar %*% sStarInv
print(variance)
print(varianceStar)

#I can calculate rough distributions of Snj given X is Unif(0,1) and gaussian kernel
#Then I can calculate an inverse given expected values
#Difficult to do all that multiplying to get variance star tho


hist(x-x0)
kappa(S[-1,-1])
kappa(sStar[-1,-1])
