SimulationData <- readRDS("Downloads/Simulation4.RData")

#Data containing errors of all iterations for the first 0th-3rd derivatives
MSEData <- SimulationData[['MSEData']]


showBoxPlots <- function(rawData, deriv = NULL,
                         interior = FALSE, showOutliers = TRUE){
  
  if(interior == TRUE){columnShift <- 4}
  else{columnShift <- 0}
  
  if(is.null(deriv) == FALSE){
    lasso <- rawData[[deriv+1]][1+columnShift]
    ridge <- rawData[[deriv+1]][2+columnShift]
    locPoly <- rawData[[deriv+1]][3+columnShift]
    
    par(mfrow=c(1,1), mar=c(0,0,0,0))
    boxplot(c(lasso, ridge, locPoly),
            main = "Comparing LPR Methods", ylab='MSE',
            names=c("Lasso", "Ridge", "LocPoly"))
  }else{
    par(mfrow=c(2,2))
    for(i in 0:3){
      nextDeriv <- c(rawData[[i+1]][1+columnShift],
                     rawData[[i+1]][2+columnShift],
                     rawData[[i+1]][3+columnShift])
      title <- paste("Derivative ", i)
      boxplot(nextDeriv, main = title, ylab='MSE', outline=showOutliers,
              names = c("Lasso", "Ridge", "LocPoly"))
    }
    par(mfrow=c(1,1))
  }
  
}

showBoxPlots(listOfSimulations[[5]][['MSEData']], NULL, TRUE, FALSE)


saveRDS(listOfSimulations, file="CompleteLassoLPRData.R")




################Violin Plots
library(ggplot2)
data <- data.frame(listOfSimulations[[4]][["MSEData"]][[1]])
data
ggplot(data, aes(x=Seed, y=LassoMSE)) + geom_violin()


















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


#Finding Min Index
which.min(unlist(listOfSimulations[[4]][["MSEData"]][[1]]['LassoMSE']))

#Finding seed given min Index
listOfSimulations[[4]][["MSEData"]][[1]]["Seed"][343,]

#Finding Error given index


