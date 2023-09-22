SimulationData <- readRDS("Downloads/Simulation4.RData")

#Data containing errors of all iterations for the first 0th-3rd derivatives
MSEData <- SimulationData[['MSEData']]


###########Determining the percentage of time Lasso was performed the best###############
HowOftenLassoBest <- rep(NULL, 4)
InteriorHowOftenLassoBest <- rep(NULL, 4)
for(j in 1:4){
  HowOftenLassoBest[j] <- prop.table(table(MSEData[[j]][,4]))[1]
  InteriorHowOftenLassoBest[j] <- prop.table(table(MSEData[[j]][,8]))[1]
}

HowOftenLassoBest
InteriorHowOftenLassoBest


showBoxPlots <- function(rawData, deriv = 0, interior = FALSE){
  if(interior == TRUE){columnShift <- 4}
  else{columnShift <- 0}
  
  allData <- unlist(c(rawData[[deriv+1]][1+columnShift],
               rawData[[deriv+1]][2+columnShift],
               rawData[[deriv+1]][3+columnShift]))
  ylim <- c(0, max(allData))
  
  par(mfrow=c(1,3))
  boxplot(rawData[[deriv+1]][1+columnShift], main='Lasso', ylim=ylim)
  boxplot(rawData[[deriv+1]][2+columnShift], main="Ridge", ylim=ylim)
  boxplot(rawData[[deriv+1]][3+columnShift], main="LocPoly", ylim=ylim)
  par(mfrow=c(1,1))
}
showBoxPlots(MSEData, 0, FALSE)

showBoxPlots(listOfSimulations[[3]][['MSEData']], 0, FALSE)

###Calculate "Average Curve" after 1000 sims

averagePoints <- rep(0, 401)

for(j in 1:401){
  allDataOfSpecPoint<- rep(0, 1000)
  for(i in 1:1000){
    allDataOfSpecPoint[i] <- SimulationData$LassoData[[i]][[1]]$LassoFit$y[j]
  }
  averagePoints[j] <- mean(allDataOfSpecPoint)
}

x <- SimulationData$LassoData[[i]][[1]]$LassoFit$x
plot(x, averagePoints, type='l')
points(x, eval(peak))
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))

SimulationData$LassoData[[1000]][[1]]$LassoFit$x[i]


######Finding Worst and Best Fit
which.max(unlist(MSEData[[1]]['LassoMSE']))
MSEData[[1]]["Seed"][88,]
#Worst fit was with seed 272955
MSEData[[1]]["LassoMSE"][88,]

#Best fit was with seed 577251
which.min(unlist(MSEData[[1]]['LassoMSE']))
MSEData[[1]]["Seed"][343,]









