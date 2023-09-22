library(glue)
RunBetterSimulation <- function(loops, equation, sampleSize, sigma,
                             seeds = as.integer(runif(loops,min=0, max=1000000))){
  #This function runs the LPRSmoothing function many times and
  #keeps track of the errors of each run. A data frame is created for the
  #original regression as well as for the 1st,2nd, and 3rd derivatives.
  #The data frames are returned in list along side a summary depicting 
  #the characteristics of this simulation
  
  columnNames <- c("LassoMSE", "RidgeMSE", "NWMSE", "LassoRank",
                   "LassoIntMSE", "RidgeIntMSE", "NWIntMSE", "LassoIntRank", "Seed")
  
  #Creates 4 empty data frames 
  ErrorDfList <- vector(mode='list', length=5)
  for(i in 1:4){
    ErrorDfList[[i]] <- data.frame(matrix(nrow=0, ncol=length(columnNames)))
    colnames(ErrorDfList[[i]]) <- columnNames
  }
  
  XData <- vector(mode='list', length=loops)
  YData <- vector(mode='list', length=loops)
  LassoOutputs <- vector(mode='list', length=loops)
  RidgeOutputs <- vector(mode='list', length=loops)
  LocPolyOutputs <- vector(mode='list', length=loops)
  
 
  
  for(i in 1:loops){
    
    errorOccured <- FALSE
    #If an error occurs due to too few datapoints within the bandwidth,
    #just return 0's for everything. will drop from dataframe post-Simulation
    out <- tryCatch( 
      {
        nextRun <- LPRSmoothingEpan(equation = equation, left = 0,
                                    right = 1, numTerms = 10,
                                    sigma = sigma, numPoints = sampleSize,
                                    seed=seeds[i])
      },
      error = function(e){
        errorOccured <<- TRUE
        print("error")
      }
    )
    
    nextRow <- vector(mode='list', length=4) #Next row to append to dataframe
    singleLassoOutput <- vector(mode='list', length=4)
    singleRidgeOutput <- vector(mode='list', length=4)
    singleLocPolyOutput <- vector(mode='list', length=4)
    
    #Error Handling
    if(errorOccured == TRUE){
      print("Error Saved")
      for(j in 1:4){
        nextRow[[j]] <- rep(0, length(columnNames))
        ErrorDfList[[j]] <- rbind(ErrorDfList[[j]], nextRow[[j]])
        colnames(ErrorDfList[[j]]) <- columnNames
      }
    }
    else{
      #Extracts data for each iteration to append into our dataframe
      for(j in 1:4){
        LassoRank <- Rank(nextRun[[j]]$LassoError, nextRun[[j]]$RidgeError,
                          nextRun[[j]]$LocpolyError)
        
        LassoIntRank <- Rank(nextRun[[j]]$LassoIntError, nextRun[[j]]$RidgeIntError,
                             nextRun[[j]]$LocpolyIntError)
        #Next row for each database
        nextRow[[j]] <- c(nextRun[[j]]$LassoError, nextRun[[j]]$RidgeError,
                          nextRun[[j]]$LocpolyError, LassoRank,
                          nextRun[[j]]$LassoIntError, nextRun[[j]]$RidgeIntError,
                          nextRun[[j]]$LocpolyIntError, LassoIntRank, seeds[i])
        
        #Appends the next row for each database
        ErrorDfList[[j]] <- rbind(ErrorDfList[[j]], nextRow[[j]])
        colnames(ErrorDfList[[j]]) <- columnNames
        
        singleLassoOutput[[j]] <- nextRun[[j]]["LassoFit"]
        singleRidgeOutput[[j]] <- nextRun[[j]]["RidgeFit"]
        singleLocPolyOutput[[j]] <- nextRun[[j]]["LocpolyFit"]
      }
      LassoOutputs[[i]] <- singleLassoOutput
      RidgeOutputs[[i]] <- singleRidgeOutput
      LocPolyOutputs[[i]] <- singleLocPolyOutput
      XData[[i]] <- nextRun[[1]]['x']
      YData[[i]] <- nextRun[[1]]['y']
    }
    #print(paste(i,"/",loops, "Loops Completed"))
    #cat("\n")
  }
  
  #If statements help judge which equation is being used to be stored in summary statement
  x <- 0.5
  if(eval(equation) == sin(5*pi*0.5)){equationName <- 'sine'}
  if(eval(equation) == 2-5*0.5 +5*exp(-400*(0.5-0.5)**2)){equationName <- 'peak'}
  if(eval(equation) == 0.3*exp(-4*(4*0.5-1)**2)+0.7*exp(-16*(4*0.5-3)**2)){equationName <- 'bimodal'}
  summary <- paste("Equation: ", equationName, " Sigma: ", sigma, " Sample Size: ", sampleSize)
  ErrorDfList[[5]] <- summary
  
  AllData <- list("MSEData" = ErrorDfList,
                  "LassoData" = LassoOutputs,
                  "RidgeData" = RidgeOutputs,
                  "LocPolyData" = LocPolyOutputs,
                  "XData" = XData,
                  "YData" = YData)
  
  return(AllData)
}

#Helper method that will compute the ranking of the Lasso compared to the other estimators
Rank <- function(lasso,ridge,locpoly){
  rank <- 1
  if(lasso > ridge){rank <- rank + 1 }
  if(lasso > locpoly){rank <- rank + 1 }
  return(rank)
}

bimodal <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2))
bugFixing <- RunBetterSimulation(3, bimodal, 100, 0.1)
bugFixing
