library('dplyr')


sampleSizes <- c(100,250,500,1000)
equations <- c("peak", "sine", "bimodal")

columns <- c("Equation", "SampleSize", "Derivative", "Seed",
             "LassoMSE", "RidgeMSE", "NWMSE",
             "LassoIntMSE", "RidgeIntMSE", "NWIntMSE",
             "LassoRank", "LassoIntRank")
combinedData <- data.frame(matrix(nrow=0, ncol=length(columns)))
combinedData
colnames(combinedData) <- columns 

for(i in 1:3){# For each Equation
  
  for(j in 1:4){ # For Each Sample Size
    
    for(k in 1:4){ #For Each Derivative
      nextDf <- data.frame(listOfSimulations[[4*(i-1)+j]]$MSEData[[k]])
      nextDf$Derivative = rep(k-1, 1000)
      nextDf$Equation = equations[i]
      nextDf$SampleSize <- sampleSizes[j]
        
      nextDf <- nextDf %>% select(Equation, SampleSize, Derivative, Seed,
                        LassoMSE, RidgeMSE, NWMSE,
                        LassoIntMSE, RidgeIntMSE, NWIntMSE,
                        LassoRank, LassoIntRank)
      
      combinedData <- rbind(combinedData, nextDf)
    }
    
  }
}
colnames(combinedData)[7] <- "LocPolyMSE"
colnames(combinedData)[10] <- "LocPolyIntMSE"
combinedData[47000:47100,]

saveRDS(combinedData, file="LPR-MSETable.RData")




combinedData$Equation <- factor(combinedData$Equation, levels=equations)
combinedData$SampleSize <- factor(combinedData$SampleSize, levels=sampleSizes)
subsection <- combinedData[combinedData$Equation =='bimodal' & combinedData$Derivative == 0, ]

ggplot(subsection, aes(x=SampleSize, y=RidgeMSE)) + geom_violin() + coord_cartesian(ylim=c(0,0.003))




