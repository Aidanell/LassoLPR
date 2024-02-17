plotBiasAndVar = function(smoothData, Data, deriv=0){
  #Regular Lasso
  lassoDf <- NULL
  for(i in 1:numExperiments){
    lassoDf <- rbind(lassoDf, Data[[i]][[deriv+1]]$LassoFit$y)
  }
  averageLassoX <- Data[[1]][[deriv+1]]$LassoFit$x #Does not change between Sims
  averageLassoY <- colMeans(lassoDf)
  
  #Smooth Lasso
  smoothLassoDf <- NULL
  for(i in 1:numExperiments){
    smoothLassoDf <- rbind(smoothLassoDf, smoothData[[i]][[deriv+1]]$LassoFit$y)
  }
  averageSmoothLassoX <- smoothData[[1]][[deriv+1]]$LassoFit$x #Does not change between Sims
  averageSmoothLassoY <- colMeans(smoothLassoDf)
  
  #Local Linear
  LLDf <- NULL
  for(i in 1:numExperiments){
    LLDf <- rbind(LLDf, Data[[i]][[deriv+1]]$LocpolyFit$y)
  }
  averageLLX <- Data[[i]][[deriv+1]]$LocpolyFit$x
  averageLLY <- colMeans(LLDf)
  
  #True Function
  x <- Data[[1]][[deriv+1]]$EvalPoints
  y<- eval(Data[[1]][[deriv+1]]$DerivFunc)
  
  par(mfrow=c(2,2))
  plot(x, y, type='l', lwd=3, main="Average Output")
  lines(averageLassoX, averageLassoY, col='red', lwd=2)
  lines(averageSmoothLassoX, averageSmoothLassoY, col='blue', lwd=2)
  lines(averageLLX, averageLLY, col='green', lwd=2)
  
  
  plot(x, abs(y - averageLassoY), type='l', col='red', main="Bias of Each Model")
  lines(x, abs(y-averageSmoothLassoY), col='blue')
  lines(x, abs(y-averageLLY), col='green')
  
  plot(x, colMeans(sweep(lassoDf, 2, averageLassoY, "-")**2), type='l', col='red', main="Variance of Each Model")
  lines(x, colMeans(sweep(smoothLassoDf, 2, averageSmoothLassoY, "-")**2), col='blue')
  lines(x, colMeans(sweep(LLDf, 2, averageLLY, "-")**2), col='green')
}
