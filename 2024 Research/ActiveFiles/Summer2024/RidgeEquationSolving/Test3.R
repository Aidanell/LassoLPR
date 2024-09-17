
n <- 20
x <- runif(n, min=-2, max=2)
sigma <- 0.05

B0 <- 3
B1 <- 0.05
y <- B0 + B1*x + rnorm(n,mean=0,sd=sigma)

X <- matrix(c(rep(1,n), x), nrow=n)
D <- matrix(c(0,0,0,1), nrow=2)
OLS <- lm(y~x)
betaOLS <- OLS$coefficients
sigmaHat <- sqrt(sum(OLS$residuals^2))

lambdaRidgeCurveB0 <- function(lambda){
  biasBetaHatRidge <- (solve(t(X) %*% X + lambda*D) %*% (t(X)%*%X %*% betaOLS))[1,1] - betaOLS[1]
  varBetaHatRidge <- sigmaHat^2 * (solve(t(X) %*% X + lambda*D) %*% (t(X)%*%X) %*% solve(t(X) %*% X + lambda*D))[1,1]
  
  MSE <- biasBetaHatRidge^2 + varBetaHatRidge
  return(MSE)
}

lambdaRidgeCurveB1 <- function(lambda){
  biasBetaHatRidge <- (solve(t(X) %*% X + lambda*D) %*% (t(X)%*%X %*% betaOLS))[2,1] - betaOLS[2]
  varBetaHatRidge <- sigmaHat^2 * (solve(t(X) %*% X + lambda*D) %*% (t(X)%*%X) %*% solve(t(X) %*% X + lambda*D))[2,2]
  MSE <- biasBetaHatRidge^2 + varBetaHatRidge
  return(MSE)
}

skip <- 0.025
input <- seq(0,25, skip)
output <- lapply(input, lambdaRidgeCurveB0)
minLambda <- skip*(which.min(output)-1)
plot(input, output)
BetaHatRidge <- solve(t(X) %*% X + minLambda*D) %*% (t(X)%*%y)
print(minLambda)

evalPoints <- seq(-2, 2, 0.1)
plot(x,y)
lines(evalPoints, BetaHatRidge[1,1] + evalPoints*BetaHatRidge[2,1])
lines(evalPoints, betaOLS[1] + evalPoints*betaOLS[2],col='red')




