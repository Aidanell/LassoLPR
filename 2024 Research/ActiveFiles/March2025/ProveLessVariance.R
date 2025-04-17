

n <- 100
x <- runif(n, min=0,1)
#x <- rnorm(n, 0.5, 0.3)
x0 = 0.5
bandwidth = 0.5

gaussKernel <- function(x){
  1/sqrt(2*pi) * exp(-x^2 / 2)
}

buildFeature <- function(midpoint, p, xValues){
  # Every j'th column is the coefficient of the (j-1)'th derivative in the
  # Taylor polynomial formed at the i'th point.
  
  X <- matrix(nrow = length(xValues), ncol = p)
  for(i in 1:p){
    X[,i] <- ((xValues - midpoint)^(i))/factorial(i)
  }
  return(X)
}

computeWeights <- function(x, x0, bandwidth, kernel='epak'){
  weights <- numeric(length(x))
  
  for(i in 1:length(x)){
    diff <- abs(x[i]-x0) / bandwidth
    if(kernel == 'epak'){
      weights[i] <- epan(diff)/bandwidth
    }else{
      weights[i] <- dnorm(diff)/bandwidth
    }
  }
  return(weights)
}


derivCalc <- function(func, numDeriv){
  drvList <- c(func)
  nextDrv <- func
  for(i in 1:numDeriv){
    nextDrv <- D(nextDrv, 'x')
    drvList <- c(drvList, nextDrv)
  }
  return(drvList)
}
derivList <- derivCalc(func = sine, numDeriv = 10)


p <- 6
termsToRemove <- c(3,4)
#For Calculating Variance
X <- buildFeature(x0, p, x)
X <- cbind(rep(1,n), X)
XStar <- X[, -(termsToRemove+1)]
W <- diag(computeWeights(x, x0, bandwidth, kernel='norm'))


sine <- expression(0.5 + sin(5*x)/3)
y <- as.matrix(eval(sine) + rnorm(n, 0, 0.1))


betahat <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y

betahatStar <- solve(t(XStar) %*% W %*% XStar) %*% t(XStar) %*% W %*% y

B0 <- eval(derivList[[1]], list(x=x0))
B1 <- eval(derivList[[2]], list(x=x0))
B2 <- eval(derivList[[3]], list(x=x0))



#Try polynomial where 3rd derivative is close to zero at x=0
testpoly <- expression(0.3*x**5 - 1.2*x**4 + 0.01*x**3 + 0.2*x**2 + 0.6*x + 0.2)

n <- 100
x0 = 0
bandwidth = 0.3
iterations <- 1000

betas <- matrix(nrow=iterations, ncol=p+1)
betaStars <- matrix(nrow=iterations, ncol=p+1)
termsToRemove <- c(3)
for(i in 1:iterations){
  x <- runif(n, -1,1)
  #x <- rnorm(n, 0.5, 0.3)
  
  X <- buildFeature(x0, p, x)
  X <- cbind(rep(1,n), X)
  XStar <- X[, -(termsToRemove+1)]
  W <- diag(computeWeights(x, x0, bandwidth, kernel='norm'))
  
  
  #sine <- expression(0.5 + sin(5*x)/3)
  y <- as.matrix(eval(testpoly) + rnorm(n, 0, 0.1))
  
  
  betahat <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
  betahatStar <- solve(t(XStar) %*% W %*% XStar) %*% t(XStar) %*% W %*% y
  
  for(j in 1:length(termsToRemove)){
    currentTerm <- termsToRemove[j]
    betahatStar <- append(betahatStar, 0, after=(currentTerm))
  }
  
  betas[i,] <- betahat
  betaStars[i,] <- betahatStar
  
}


var(betas[,2])
var(betaStars[,2])

hist(betas[,2])
hist(betaStars[,2])

polyderivList <- derivCalc(func = testpoly, numDeriv = 10)
B0 <- eval(polyderivList[[1]], list(x=x0))
B1 <- eval(polyderivList[[2]], list(x=x0))

for(i in 1:10){
  print(eval(polyderivList[[i]], list(x=x0)))
}
j <- 5
beta <- eval(polyderivList[[j+1]], list(x=x0))
betahat <- mean(betas[,j+1])
betahatStar <- mean(betaStars[,j+1])
print(beta)
print(betahat)
print(betahatStar)
var(betas[,j+1])
var(betaStars[,j+1])

xlimits <- c(min(betas[,j+1]), max(betas[,j+1]))
hist(betas[,j+1], xlim=xlimits, breaks=20)
hist(betaStars[,j+1], xlim=xlimits, breaks=20)
