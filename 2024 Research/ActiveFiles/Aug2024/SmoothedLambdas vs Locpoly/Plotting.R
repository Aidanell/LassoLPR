

deriv <- 4

locpolyBandwidth1 <- thumbBw(x,y,0, gaussK)
locpolyBandwidth2 <- thumbBw(x,y,3, gaussK)

locpolyFit1 <- KernSmooth::locpoly(x,y,degree=deriv+1,drv=deriv, bandwidth=locpolyBandwidth1, gridsize=401)
locpolyFit2 <- KernSmooth::locpoly(x,y,degree=deriv+1,drv=deriv, bandwidth=locpolyBandwidth2, gridsize=401)




temp <- x
x <- evalPoints

par(mfrow=c(1,1))

plot(evalPoints, eval(derivList[deriv+1]), type='l', lwd=2)
#lines(evalPoints, smoothLassoOutput[,deriv+1], col="orange", lwd=2)
lines(locpolyFit1, lwd=2, col='blue')
lines(locpolyFit2, lwd=2, col='red')
#lines(derEst$x, derEst$der)
x <- temp


