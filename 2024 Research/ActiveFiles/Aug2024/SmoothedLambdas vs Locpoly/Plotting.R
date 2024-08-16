

deriv <- 2

ThumbBandwidth <- thumbBw(x,y,deg=deriv, gaussK)

locpolyFit <- KernSmooth::locpoly(x,y,degree=deriv+1,drv=deriv, bandwidth=ThumbBandwidth, gridsize=401)
#derEst <- compDerEst(x,y,deriv-1)



temp <- x
x <- evalPoints

par(mfrow=c(1,1))

plot(evalPoints, eval(derivList[deriv+1]), type='l', lwd=2)
lines(evalPoints, smoothLassoOutput[,deriv+1], col="orange", lwd=2)ww
lines(locpolyFit, lwd=2, col='blue')
lines(derEst$x, derEst$der)
x <- temp


