library(locpol)
library(KernSmooth)

#Create neccesary data----
set.seed(100)
n <- 500
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
#peak <- expression(0.5 + sin(5*x)/3)
#peak <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2))
p <- 9
sigma <- sqrt(0.5)
#sigma <- 0.1
deriv <- 1#Used for bandwidth calculations
sampleData <- buildData(n, peak, sigma)
bandwidth <- thumbBw(sampleData$x,sampleData$y,deriv,gaussK)


#Builds the true datapoints for the function. Helps with MSE----
derivList <- derivCalc(func = peak, numDeriv = 10)
x <- seq(0,1,length.out=401)
trueFunctions <- list(eval(derivList[1]),eval(derivList[2]),eval(derivList[3]),eval(derivList[4]), eval(derivList[5]),
                      eval(derivList[6]), eval(derivList[7]),eval(derivList[8]),eval(derivList[9]),eval(derivList[10]))
#Calculate true Y for simulated x's
x <- sampleData$x
trueY <- eval(peak)



#Running the Method Simulation----

customCVLasso <- CustomLassoLPR(sampleData$x, sampleData$y, h=bandwidth, p=p, varBands=FALSE)

#Old LassoLPR
LPRband <-thumbBw(sampleData$x,sampleData$y,3,EpaK)
lassoLPR <- lassoLPR(sampleData$x, sampleData$y, h=LPRband, p=p)


derivative <- 0
gridPoints <- seq(0, 1, length.out=401)

#Locpoly Comparison
locBand <- dpill(sampleData$x, sampleData$y)
locOut <- locpoly(sampleData$x, sampleData$y, bandwidth=locBand*3, drv=derivative)
#plotting Lambdas
plot(gridPoints, customCVLasso[[4]], type='l', col='blue')
lines(gridPoints, lassoLPR$lambdas, type='l', col='green')

#Plotting Derivative Comparisons
plot(gridPoints, trueFunctions[[derivative+1]], type='l', lwd=2, ylab='y', xlab='x')
#lines(locOut, col='purple', lwd=2) #Locpoly Example
#lines(gridPoints, customCVLasso[[1]][,derivative+1], col='blue', type='l', lwd=2) #Unsmoothed Custom CV
lines(gridPoints, lassoLPR$lasso[,derivative+1], col='green', type='l', lwd=2) #Old Smoothed LassoLPR
lines(gridPoints, customCVLasso[[3]][,derivative+1], col='red', type='l', lwd=2) #Smoothed Custom CV


#lines(gridPoints, customCVLasso[[4]], col='red')



