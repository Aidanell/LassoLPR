#############Organizing Data#################
#Getting the percentage of times Lasso was the leading method for each simulation
DataMatrixInt <- matrix(nrow=12, ncol=4)
for(i in 1:12){
  for(j in 1:4){
    DataMatrixInt[i,j] <- prop.table(table(listOfSimulations[[i]][["MSEData"]][[j]][,4]))['1']
  }
}

colnames(DataMatrixInt) <- c("0th Deriv", "1st Deriv", "2nd Deriv", "3rd Deriv")
rownames(DataMatrixInt) <- c("Peak: 100", "Peak: 250", "Peak: 500", "Peak: 1000",
                          "Sine: 100", "Sine: 250", "Sine: 500", "Sine: 1000",
                          "Bimodal: 100", "Bimodal: 250", "Bimodal: 500", "Bimodal: 1000")
DataMatrixInt
view(DataMatrix)
saveRDS(DataMatrix, file="PercentageLassoDidBestTable.R")




calcLassoPerformance <- function(simData, interior = FALSE){
  DataMatrix <- matrix(nrow=12, ncol=4)
  if(interior == TRUE){rankRow <- 8}
  else{rankRow <- 4}
  
  for(i in 1:12){
    for(j in 1:4){
      DataMatrix[i,j] <- prop.table(table(simData[[i]][["MSEData"]][[j]][,rankRow]))['1']
    }
  }
  
  colnames(DataMatrix) <- c("0th Deriv", "1st Deriv", "2nd Deriv", "3rd Deriv")
  rownames(DataMatrix) <- c("Peak: 100", "Peak: 250", "Peak: 500", "Peak: 1000",
                               "Sine: 100", "Sine: 250", "Sine: 500", "Sine: 1000",
                               "Bimodal: 100", "Bimodal: 250", "Bimodal: 500", "Bimodal: 1000")
  
  return(DataMatrix)
}

lassoPerformTable<- calcLassoPerformance(listOfSimulations, FALSE)
lassoIntPerformTable <- calcLassoPerformance(listOfSimulations, TRUE)
saveRDS(lassoPerformTable, file="lassoPerformTable.R")
saveRDS(lassoIntPerformTable, file="lassoIntPerformTable.R")







#This is taking data and comparing it to each method based on percentage difference in error
ComparingMethodList <- vector(mode='list', length=4)
for(i in 1:4){
  ComparingMethodList[[i]] <- matrix(nrow=12, ncol=3)
  colnames(ComparingMethodList[[i]]) <- c("Ridge", "Nadaraya-Watson", "Local Linear")
  rownames(ComparingMethodList[[i]]) <- c("Peak: 100", "Peak: 250", "Peak: 500", "Peak: 1000",
                                          "Sine: 100", "Sine: 250", "Sine: 500", "Sine: 1000",
                                          "Bimodal: 100", "Bimodal: 250", "Bimodal: 500", "Bimodal: 1000")
  }

for(i in 1:4){#Each Derivative
  for(j in 1:12){ #Each simulation
    for(k in 2:4){ #Each method
      ComparingMethodList[[i]][j,k-1] <- sapply((SimData[[j]][[i]][k] - SimData[[j]][[i]][1]) / SimData[[j]][[i]][1], mean, na.rm = TRUE)
    }
  }
}
ComparingMethodList

saveRDS(ComparingMethodList, file="PercentDifferenceTables.R")

