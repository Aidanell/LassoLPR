plotBiasAndVar = function(SimData, deriv=0){
  #Regular Lasso
  lassoDf <- NULL
  for(i in 1:length(SimData)){
    lassoDf <- rbind(lassoDf, SimData[[i]]$lassoOutput[,deriv+1])
  }
  
  averageLassoY <- colMeans(lassoDf)
  
  #Smooth Lasso
  smoothLassoDf <- NULL
  for(i in 1:numExperiments){
    smoothLassoDf <- rbind(smoothLassoDf, SimData[[i]]$smoothLassoOutput[,deriv+1])
  }
  averageSmoothLassoY <- colMeans(smoothLassoDf)
  
  
  #Local Linear
  LLDf <- NULL
  for(i in 1:numExperiments){
    LLDf <- rbind(LLDf, SimData[[i]]$locpolyOutput[,deriv+1])
  }
  averageLLY <- colMeans(LLDf)
  
  #True Function
  x <- SimData[[1]]$evalPoints
  y<- eval(SimData[[1]]$derivFunctions[deriv+1])
  
  if(length(y) < 2){
    y <- rep(y, 401)
  }
  
  #Bias Calculations
  lassoBias <- abs(y - averageLassoY)
  smoothLassoBias <- abs(y-averageSmoothLassoY)
  locpolyBias <- abs(y-averageLLY)
  #Variance Calculations
  lassoVar <- colMeans(sweep(lassoDf, 2, averageLassoY, "-")**2)
  smoothLassoVar <- colMeans(sweep(smoothLassoDf, 2, averageSmoothLassoY, "-")**2)
  locpolyVar <- colMeans(sweep(LLDf, 2, averageLLY, "-")**2)
  #MSE Calculations
  lassoMSE <- lassoVar + lassoBias**2
  smoothLassoMSE <- smoothLassoVar + smoothLassoBias**2
  locpolyMSE <- locpolyVar + locpolyBias**2

  #plotting
  #par(mfrow=c(2,2))
  
  layout(matrix(c(1,2,3,3,4,5), 3, 2, byrow=TRUE),
         heights=c(3,1,3))
  par(mar=c(2.5,2,2,1))
  
  combinedOutputs <- c(averageLassoY, averageSmoothLassoY, y)
  plotLimits <- c(min(combinedOutputs), max(combinedOutputs))
  plot(x, y, type='l', lwd=3, main="Average Output", ylim=plotLimits, cex.main=1.5)
  lines(x, averageLassoY, col='red', lwd=2)
  lines(x, averageSmoothLassoY, col='orange', lwd=2)
  lines(x, averageLLY, col='blue', lwd=2)
  
  #Bias
  combinedOutputsBias <- c(lassoBias, smoothLassoBias)
  plotLimitsBias <- c(min(combinedOutputsBias), max(combinedOutputsBias))
  plot(x, lassoBias, type='l', col='red', main="Bias of Each Model",
       lwd=2, ylab="Bias", ylim=plotLimitsBias, cex.main=1.5)
  lines(x, smoothLassoBias, col='orange', lwd=2)
  lines(x, locpolyBias, col='blue', lwd=2)
  
  
  #Legend
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legendTitle = paste("Derivative: ", deriv)
  legend("bottom", legend=c("True Function", "Lasso", "Locpoly", "Smoothed Lambdas Lasso"),
         col=c("black", "red", "blue", "orange"), pch=20, pt.cex=4, horiz=TRUE, cex=1.3)
  mtext(legendTitle, at=0.5, cex=1.5)
  

  #Variance
  combinedOutputsVar <- c(lassoVar, smoothLassoVar)
  plotLimitsVar <- c(min(combinedOutputsVar), max(combinedOutputsVar))
  plot(x, lassoVar, type='l', col='red', main="Variance of Each Model",
       lwd=2, ylab="Variance", ylim=plotLimitsVar, cex.main=1.5)
  lines(x, smoothLassoVar, col='orange', lwd=2)
  lines(x, locpolyVar, col='blue', lwd=2)
  
  #MSE
  plot(x, lassoMSE, type='l', col='red', main="MSE of Each Model",
       lwd=2, ylab="MSE", cex.main=1.5)
  lines(x, smoothLassoMSE, col='orange', lwd=2)
  lines(x, locpolyMSE, col='blue', lwd=2)
  
}
  
  
