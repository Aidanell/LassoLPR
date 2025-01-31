library(locpol)
library(KernSmooth)

#Create neccesary data----
set.seed(200)
n <- 500
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))  # sigma = sqrt(0.5)
#peak <- expression(0.5 + sin(5*x)/3)
#peak <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2))
p <- 13
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

#Calculate Variable Bandwidth Coefficients----
#This section roughly calculates the variability of the data around that point. The rationale is to
#reduce bandwidth in higher variability areas. I then linearly transform the variance to be in range of [0,3]
#NOTE: The transform is function specific, and I choose this as a starting point for the peak function,
#This is very far from a rule of thumb to use generally

gridPoints <- seq(0,1, length.out=401)
localVariance <- numeric(401)
windowSize <- 0.1
for(i in 1:length(gridPoints)){
  currPoint <- gridPoints[i]
  InWindow <- which( abs(sampleData$x - currPoint)/windowSize < 1)
  meanInWindow <- mean(sampleData$y[InWindow]^2)
  localVariance[i] <- mean(sampleData$y[InWindow]^2) - mean(sampleData$y[InWindow])^2
}
plot(gridPoints, localVariance) #Show Variance
variableBands <- 3 - localVariance/3 #Transform data to represent multipliers for the bandwidth
plot(gridPoints, variableBands, type='l')


#Running the Method Simulation----
customCVLasso <- testingCustomCV(sampleData$x, sampleData$y, h=bandwidth, p=p, varBands=FALSE)

LPRband <-thumbBw(sampleData$x,sampleData$y,3,EpaK)
lassoLPR <- lassoLPR(sampleData$x, sampleData$y, h=LPRband, p=p)

#Plotting----
#Choose what derivative to plot, Change as needed*********

plotOutput <- function(CVObject, lassoLPRObject, derivative, sampleData, plotLam=FALSE){
  
  if(plotLam == FALSE){
    #Runs a locpoly example to compare with.
    locBand <- dpill(sampleData$x, sampleData$y)
    locOut <- locpoly(sampleData$x, sampleData$y, bandwidth=locBand*(derivative+1), drv=derivative)
    
    
    title <- paste("Performance Comparison of Lasso vs Locpoly for Derivative ", derivative+1)
    plot(gridPoints, trueFunctions[[derivative+1]], type='l', lwd=2, ylab='y', xlab='x', main=title)
    #lines(locOut, col='red', lwd=2)
    lines(gridPoints, CVObject[[1]][,derivative+1], col='blue', type='l', lwd=2)
    lines(gridPoints, lassoLPRObject[[1]][,derivative+1], col='green', type='l', lwd=2)
    lines(gridPoints, CVObject[[3]][,derivative+1], col='red', type='l', lwd=2)
    #Legend
    #legend('topleft', legend=c('Custom CV-LassoLPR', "Old LassoLPR", "True Function"), fill=c("blue", "green", "black"))
    #Optionally plot with Datapoints (only relevant for derivative=0)
    if(derivative == 0){
      points(sampleData$x, sampleData$y, pch=20, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1))  
    }
  }else{
    plot(gridPoints, CVObject[[2]], type='l', col='blue', main="Lambdas of CustomCV vs lassoLPR")
    lines(gridPoints, lassoLPRObject$lambdas, col='green')
    #lines(gridPoints, lassoLPRObject$unsmoothLambdas, col="purple")
    lines(gridPoints, CVObject[[4]], col="red")
  }
}

plotOutput(customCVLasso, lassoLPR, 2, sampleData, plotLam=TRUE)
plotOutput(customCVLasso, lassoLPR, 2, sampleData, plotLam=FALSE)













