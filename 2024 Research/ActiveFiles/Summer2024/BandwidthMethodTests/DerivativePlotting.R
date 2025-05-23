temp <- x
x <- evalPoints
deriv <- 3
#plot(evalPoints, eval(derivList[deriv+1]), main="Lasso Spatial (Red) vs Ridge Spatial (Blue), 1st Derivative")
#lines(evalPoints, SpatialLassoOutput[,deriv+1], col='red', lwd=2)
#lines(evalPoints, SpatialRidgeOutput[,deriv+1], col='blue', lwd=2)
#lines(evalPoints, testLassoEpan$smoothLassoOutput[,deriv+1], col="orange", lwd=2)
#lines(evalPoints, testLassoEpan$lassoOutput[,deriv+1], col="purple", lwd=2)

par(mfrow=c(1,1))
plot(evalPoints, eval(derivList[deriv+1]), main="Lasso Thumb (Red) vs Ridge Thumb (Blue), 1st Derivative")
#lines(evalPoints, lassoOutput[,deriv+1], col='red', lwd=2)
#lines(evalPoints, ridgeOutput[,deriv+1], col='blue', lwd=2)
lines(evalPoints, testLassoEpan$smoothLassoOutput[,deriv+1], col="orange", lwd=2)
#lines(evalPoints, testLassoEpan$lassoOutput[,deriv+1], col="purple", lwd=2)
x <- temp


#Plot lambdas
#plot(evalPoints, lambdaTracking)
#points(evalPoints, testLassoEpan$lambdas$smoothed, col='orange')
#points(evalPoints, testLassoEpan$lambdas$unsmoothed, col='purple')



#plot(x,y)
#lines(evalPoints, lassoOutput[,deriv+1], col='red', lwd=2)
#lines(evalPoints, ridgeOutput[,deriv+1], col='blue', lwd=2)
#lines(evalPoints, testLassoEpan$smoothLassoOutput[,1], col="orange", lwd=2)
#lines(evalPoints, testLassoEpan$lassoOutput[,1], col="purple", lwd=2)
