#The Simulation runs 100 epochs for each section
#Our bandwidth is different for each derivative we want to specify
#so we test at derivatives 0,1,and 3. 
#Locpoly is done the same way and compared with
library(glmnet)
library(KernSmooth)
library(stats)
library(epandist)
library(locpol)
#Initially trying for 100 epochs at n=500
#Publishable numbers should prob have 1000 epochs and n=c(100,500,1000)

runSection <- function(epochs, n,p,equation, sigma, deriv){
  MasterDataList <- vector("list", epochs)
  for(i in 1:epochs){
    MasterDataList[[i]] <- lassoEpoch(n, p, equation, sigma, deriv)
  }
  return(MasterDataList)
}

#Universal Parameters
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2)) # sigma = sqrt(0.5)
equation <- peak
n <- 500
p <- 10
sigma <- sqrt(0.5)
gridPoints <- seq(0,1,length.out=401)
epochs <- 250

LassoD0 <- runSection(epochs, n, p, equation, sigma, deriv=0)
LocoutD0 <- locpolyEpochs(epochs, n, equation, sigma, deriv=0)
print("D0 Done!")
LassoD1 <- runSection(epochs, n, p, equation, sigma, deriv=1)
LocoutD1 <- locpolyEpochs(epochs, n, equation, sigma, deriv=1)
print("D1 Done!")
LassoD3 <- runSection(epochs, n, p, equation, sigma, deriv=3)
LocoutD3 <- locpolyEpochs(epochs, n, equation, sigma, deriv=3)
print("D3 Done!")

#For Plotting Purposes Only-----
aveLassoD0 <- calcAverage(LassoD0, deriv=0)
aveLocOutD0 <- colMeans(LocoutD0)
aveLassoD1 <- calcAverage(LassoD1, deriv=1)
aveLocOutD1 <- colMeans(LocoutD1)
aveLassoD3 <- calcAverage(LassoD3, deriv=3)
aveLocOutD3 <- colMeans(LocoutD3)

#True Output
derivList <- derivCalc(func = equation, numDeriv = 5)
x <- gridPoints
trueFunctions <- list(eval(derivList[1]),eval(derivList[2]),eval(derivList[3]),eval(derivList[4]), eval(derivList[5]),eval(derivList[6]))

LassoVarD0 <- calcVariance(LassoD0, deriv=0)
LocVarD0 <- calcLocpolyVar(LocoutD0, trueFunctions[[1]])
LassoVarD1 <- calcVariance(LassoD1, deriv=1)
LocVarD1 <- calcLocpolyVar(LocoutD1, trueFunctions[[2]])
LassoVarD3 <- calcVariance(LassoD3, deriv=3)
LocVarD3 <- calcLocpolyVar(LocoutD3, trueFunctions[[4]])

#Plotting True Function
plot(gridPoints, trueFunctions[[1]], main="True Function vs Average Output", ylab='y', xlab=" ", type='l', lwd=2)
lines(gridPoints, aveLassoD0$aveUnsmoothOutput, col='red')
lines(gridPoints, aveLassoD0$aveSmoothOutput, col='orange')
lines(gridPoints, aveLocOutD0, col='blue')

#Plotting Variance
plot(gridPoints, LassoVarD0$varUnsmoothOutput, type='l', col='red', main="Variance of Methods", ylab='Variance', xlab=" ")
lines(gridPoints, LassoVarD0$varSmoothOutput, col='orange')
lines(gridPoints, LocVarD0, col='blue')

#Plotting Bias
unsmoothLassoBias <- abs(aveLassoD0$aveUnsmoothOutput - trueFunctions[[1]])
smoothLassoBias <- abs(aveLassoD0$aveSmoothOutput - trueFunctions[[1]])
locOutBias <- abs(aveLocOutD0 - trueFunctions[[1]])
plot(gridPoints, unsmoothLassoBias, col='red', type='l', main="Bias of Methods", ylab="Bias", xlab=" ")
lines(gridPoints, smoothLassoBias, col='orange')
lines(gridPoints, locOutBias, col='blue')

#Plotting MSE
plot(gridPoints, unsmoothLassoBias**2 + LassoVarD0$varUnsmoothOutput, col='red', type='l', main="MSE of Methods", ylab="MSE", xlab=" ")
lines(gridPoints, smoothLassoBias**2 + LassoVarD0$varSmoothOutput, col='orange')
lines(gridPoints, locOutBias**2 + LocVarD0, col='blue')

mtext("Simulation Results for 0th Derivative, 250 epochs, n=500", side = 3, line = -21, outer = TRUE)
mtext("LassoUnsmooth = Red, LassoSmooth = Orange, Locpoly = Blue", side = 3, line = -22, outer = TRUE)






