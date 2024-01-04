FindOptimalBand <- function(
    equation,
    sigma,
    sampleSize,
    seed = NULL
    ){
  bandwidthMultipliers <- seq(1,15,0.25) #Bandwidths to evaluate LPR
  if(is.null(seed) == TRUE){
    seed = runif(1, min=0, max=1000)
  }
  LassoErr0 <- LassoErr1 <- LassoErr2 <- LassoErr3 <- c()
  for(i in 1:length(bandwidthMultipliers)){
    
    nextRun <- LPRSmoothing(
      equation = equation,
      left = 0,
      right = 1,
      numTerms = 10,
      sigma = sigma,
      numPoints = sampleSize,
      plot = TRUE,
      seed = seed,
      bandMulti = bandwidthMultipliers[i]
    )
    
    LassoErr0 <- c(LassoErr0, nextRun[[1]]$LassoError)
    LassoErr1 <- c(LassoErr1, nextRun[[2]]$LassoError)
    LassoErr2 <- c(LassoErr2, nextRun[[3]]$LassoError)
    LassoErr3 <- c(LassoErr3, nextRun[[4]]$LassoError)
  }
  
  dev.off()
  par(mfrow=c(2,2))
  plot(bandwidthMultipliers, LassoErr0, main="Lasso-LPR ")
  plot(bandwidthMultipliers, LassoErr1)
  plot(bandwidthMultipliers, LassoErr2)
  plot(bandwidthMultipliers, LassoErr3)
  
  df <- data.frame(Deriv0 = LassoErr0, Deriv1 = LassoErr1,
                   Deriv2 = LassoErr2, Deriv3 = LassoErr3,
                   bandwidth = bandwidthMultipliers)
  print(df)
  return(df)
}




MultipleRuns <- function(iterations, equation, sigma, sampleSize){
  
  MSEData <- vector(mode='list', length=iterations)
  for(i in 1:iterations){
    nextDf <- FindOptimalBand(equation, sigma, sampleSize)
    MSEData[[i]] <- nextDf
    print(paste(i, " Loop Done"))
  }
  par(mfrow=c(1,1))
  
  #Plotting
  plot(MSEData[[1]]$bandwidth, MSEData[[1]]$Deriv0, ylim=c(0,max(MSEData[[1]]$Deriv0)))
  for(i in 2:iterations){
    lines(MSEData[[i]]$bandwidth, MSEData[[i]]$Deriv0)
  }
  return(MSEData)
}


peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))
sine <- expression(sin(5*pi*x))
bimodal <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2))


plot(BandTest3[[1]]$bandwidth, BandTest3[[1]]$Deriv0, type='l', ylim=c(0,0.001))
for(i in 1:5){
  lines(BandTest3[[i]]$bandwidth, BandTest3[[i]]$Deriv0)
}
BandTest2 <- MultipleRuns(5, peak, sqrt(0.5), 500)
print("Test 2 Complete")
BandTest3 <- MultipleRuns(5, bimodal, 0.1, 500)
print("Test 3 Complete")
BandTest4 <- MultipleRuns(5, peak, sqrt(0.5), 100)
print("Test 4 Complete")
BandTest5 <- MultipleRuns(5, peak, sqrt(0.5), 1000)
print("Test 5 Complete")




#All below was me finding the average bandwidth that minimized the MSE.
#Final result over 25 iterations with varying equations and sample size was
#8.22

print(5 + 0.5 * which.min(BandTest2[[1:5]]$Deriv0))

Test1Deriv0 <- c()
for(i in 1:5){
  Test1Deriv0 <- c(Test1Deriv0, 5 + 0.5 * which.min(BandTest1[[i]]$Deriv0))
}
Test1Deriv0

Test2Deriv0 <- c()
for(i in 1:5){
  Test2Deriv0 <- c(Test2Deriv0, 5 + 0.5 * which.min(BandTest2[[i]]$Deriv0))
}
Test2Deriv0

Test3Deriv0 <- c()
for(i in 1:5){
  Test3Deriv0 <- c(Test3Deriv0, 5 + 0.5 * which.min(BandTest3[[i]]$Deriv0))
}
Test3Deriv0

Test4Deriv0 <- c()
for(i in 1:5){
  Test4Deriv0 <- c(Test4Deriv0, 5 + 0.5 * which.min(BandTest4[[i]]$Deriv0))
}
Test4Deriv0

Test5Deriv0 <- c()
for(i in 1:5){
  Test5Deriv0 <- c(Test5Deriv0, 5 + 0.5 * which.min(BandTest5[[i]]$Deriv0))
}
Test1Deriv0

TestMeans <- c(mean(Test1Deriv0), mean(Test2Deriv0), mean(Test3Deriv0), mean(Test4Deriv0), mean(Test5Deriv0))
TestMeans
mean(TestMeans)
