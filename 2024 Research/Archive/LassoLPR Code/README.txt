Here is a quick rundown about my code.

The "1.LPRSmoothing.R" is the most important part. It runs a single iteration of creating data, applying Lasso-LPR and other non-parametric regression methods.
It then plots the data cleanly and returns all the data in a neat format. This file is by far the most important to be reviewed. There is also a few helper methods which are called repeatedly in the LPRSmoothing method.

"2.RunOneSimulation.R" wraps the LPRSmoothing method and will repeatedly simulate the method with the same parameters but different data generated. 
It returns a data frame where each row is an iteration and columns include the error for each method and where the Lasso ranks out of the 4 methods used (Lasso-LPR, Ridge-LPR, Nadaraya-Watson, Local Linear).
A calculation of the error when factoring out the edges of the data for each method is also included. This also contains some error catching in the off chance no data points are found within the bandwidth.
In my large simulation of 12,000 iterations, this only occured twice but had to accounted for otherwise the entire simulation would fail.

"3.AllSimulations.R" wraps RunOneSimulation and will run that method for every combination of regression function and sample size desired. In total, 12,000 iterations of the LPRSmoothing method will be used so I do not reccomend
running this function. This function is the simulation I ran which took ~6 days to complete. It also returns the data in a fairly neat way.

"4.SimDataAnalysis.R" This is a simple file I used to clean up the data from AllSimulations.R and output a few cleaner tables that are easy to read

"5.FindOptimalBand.R" This was my simulation to find the bandwidth multiplier that worked best. It will take the same data and evaluate it at many different bandwidths and return a list of the MSE's at each bandwidth.
The MultipleRuns functions will run different seeds with the same equation, sigma, and sample size, which can all be plotted easily. This is also takes a while to run so I do not recommend executing this file.


"LPRSimulationData" is just the raw data returned from running AllSimulations.R

"PercentageLassoDidBestTable" tells you how often Lasso did best for each scenario. For example, a value of 0.949 means in 949/1000 iterations, Lasso-LPR had lower MSE than ALL of the other 3 methods used.

"PercentageDifferenceTables" tells you the average percentage difference the MSE was for each method. There a four tables which are in order of the original function, then its first 3 derivatives.
For example, a value of 1.84 means Lasso-LPR on average had 184% less error than the method specified in the column


Other Notes:
-The "1.LPRSmoothing.R" is really the only code that needs reviewing. The others are primarily just wrappers that helped me run iterations many times and organize the data easier.

Please email me if you have any other questions about the code.
Thanks,
-Aidan