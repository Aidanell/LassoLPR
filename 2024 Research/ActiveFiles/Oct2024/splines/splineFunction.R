
library(locpol)
library(lassoLPR)
buildData <- function(n, equation, sigma){
  x <- sort(runif(n, min=0, max=1))
  noise <- rnorm(n = length(x), mean = 0, sd = sigma)
  trueY <- eval(equation)
  y <- trueY + noise
  return(list(x=x,y=y,trueY=trueY))
}


findKnots <- function(lassoData, deriv=3, quant = 0.5){
  critIndicies <- which(abs(diff(sign(diff(lassoData$lasso[,deriv+1]))))==2)+1
  critPoints <- gridPoints[critIndicies]
  
  #Use Dilation to check where local mins and maxes are
  maxDilation <- VPdtw::dilation(lassoData$lasso[,deriv+1], 5)
  minDilation <- VPdtw::dilation(-lassoData$lasso[,deriv+1], 5)
  dilationInvariant <- which(lassoData$lasso[,deriv+1] == maxDilation | -lassoData$lasso[,deriv+1] == minDilation)
  
  #Only points that are critical and dilation invariant
  knotIndex <- critIndicies[which(critIndicies %in% dilationInvariant)]
  knots <- gridPoints[knotIndex]
  knotSize <- lassoData$lasso[knotIndex,deriv+1]
  
  #Get rid of the lowest percentage knots in magnitude
  threshhold <- quantile(abs(knotSize), quant)
  knotIndex <- knotIndex[which(abs(knotSize) > threshhold)]
  knots <- gridPoints[knotIndex]
  
  return(list("knots"=knots, "knotIndex"=knotIndex))
}
#Create neccesary data----
n <- 250
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
#peak <- expression(0.5*x)
p <- 9
sigma <- sqrt(0.5)
deriv <- 3 #Used for bandwidth calculations
sampleData <- buildData(n, peak, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv,EpaK)
gridPoints <- seq(0,1,length.out=401)
x <- seq(0,1,length.out=401)

testLasso <- lassoLPR(sampleData$x, sampleData$y, bandwidth,p)

knotData3 <- findKnots(testLasso, 3, 0.6)
knotData2 <- findKnots(testLasso, 2, 0.6)

knotIndex <- sort(c(knotData2$knotIndex, knotData3$knotIndex))
knots <- sort(c(knotData2$knots, knotData3$knots))

#Getting Spline of 2nd/3rd deriv knots
dataFrame <- data.frame(sampleData)
splineOutput <- lm(y ~ splines::ns(x, knots=knots), data=dataFrame)
splinePred <- predict(splineOutput, data.frame(gridPoints))

splineOutput2 <- lm(y ~ splines::ns(x, knots=knotData2$knots), data=dataFrame)
splinePred2 <- predict(splineOutput2, data.frame(gridPoints))

splineOutput3 <- lm(y ~ splines::ns(x, knots=knotData3$knots), data=dataFrame)
splinePred3 <- predict(splineOutput3, data.frame(gridPoints))

plot(sampleData, main="Spline with knots at 2nd Derivative Maximums")
lines(gridPoints, splinePred2, col='purple', lwd=2)

plot(sampleData, main="Spline with knots at 3rd Derivative Maximums")
lines(gridPoints, splinePred3, col='purple', lwd=2)

plot(sampleData, cex=0.25, main="Spline with knots at 2nd and 3rd Derivative Maximums")
lines(gridPoints, splinePred, col='purple', lwd=2)
points(knots, splinePred[knotIndex], col='red' , pch=16)
lines(gridPoints, testLasso$lasso[,1], col='blue' ,lwd=2)
lines(gridPoints, trueFunctions[[1]], lwd=2)



plot(gridPoints, testLasso$lasso[,3], type='l')
points(knotData2$knots, testLasso$lasso[knotData2$knotIndex,3], col='red')



#Use Box-Cox
#Use Q1y and Q2y with QQ-plot