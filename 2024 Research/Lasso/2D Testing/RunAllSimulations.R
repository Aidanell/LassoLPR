
#This function runs the RunOneSimulation function repeatedly for
#every possible scenario we want to test. It returns a list of returned lists
#from the RunOneSimulation function.
iterations <- 1000
  
#Desired equations
peak <- expression(2-5*x +5*exp(-400*(x-0.5)**2))
sine <- expression(sin(5*pi*x))
bimodal <- expression(0.3*exp(-4*(4*x-1)**2)+0.7*exp(-16*(4*x-3)**2))
  
equations <- c(peak, sine, bimodal)
sigmas <- c(sqrt(0.5), 0.5, 0.1) #Sigmas that are appropriate for their respective functions
sampleSizes <- c(100, 250, 500, 1000) #Each function will be tested for each sample size
  
seeds <- as.integer(runif(iterations,min=0, max=1000000))
  
listOfSimulations <- vector(mode="list", length=12) #Each entry contains data from each sim
for(i in 1:3){
  currentEquation <- equations[i]
  currentSigma <- sigmas[i]
  for(j in 1:4){
    currentSampleSize <- sampleSizes[j]
    nextSim <- RunBetterSimulation(
                                loops = iterations,
                                equation = currentEquation,
                                sampleSize = currentSampleSize,
                                sigma = currentSigma,
                                seeds = seeds)
    
    filename <- glue("researchData/Simulation{4*(i-1)+j}.RData")
    saveRDS(nextSim, file=filename)
    listOfSimulations[[4*(i-1)+j]] <- nextSim #Append Simulation to correct spot
    print(paste("Sim ", 4*(i-1)+j, "Complete"))
    
  }
}
  
listOfSimulations




