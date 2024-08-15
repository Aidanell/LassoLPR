x <- c(0,1,2)
y <- c(4,5,6)
n <- length(x)

X <- matrix(c(rep(1,n), x), nrow=n)
D <- matrix(c(0,0,0,1), nrow=2)
TrueBeta <- matrix(c(4,1), nrow=2)
lambda <- 0.5
#BetaRidge = (T(X)X + lambda*D)^-1 * T(X)*y


BetaHatRidge <- solve(t(X) %*% X + lambda*D) %*% t(X)%*%y
BetaHat <- solve(t(X) %*% X) %*% t(X)%*%y
plot(x,y)
abline(BetaHatRidge)
abline(BetaHat, col='blue')

EbetaHatRidge <- solve(t(X) %*% X + lambda*D) %*% (t(X)%*%X %*% TrueBeta)
VbetaHatRidge <- solve(t(X) %*% X + lambda*D) %*% (t(X)%*%X) %*% solve(t(X) %*% X + lambda*D)

EbetaHat <- TrueBeta
VbetaHat <- solve(t(X) %*% X)

curve(dnorm(x, EbetaHatRidge[1,1], VbetaHatRidge[1,1]), from=1, to=7)
curve(dnorm(x, EbetaHat[1,1], VbetaHat[1,1]), from=1, to=7, add=TRUE, col="blue")

curve(dnorm(x, EbetaHatRidge[2,1], VbetaHatRidge[2,2]), from=-1, to=2)
curve(dnorm(x, EbetaHat[2,1], VbetaHat[2,2]), from=-1, to=2, add=TRUE, col="blue")

pnorm(1.5, EbetaHatRidge[2,1], VbetaHatRidge[2,2]) - pnorm(0.5, EbetaHatRidge[2,1], VbetaHatRidge[2,2])
pnorm(1.5, EbetaHat[2,1], VbetaHat[2,2]) - pnorm(0.5, EbetaHat[2,1], VbetaHat[2,2])
