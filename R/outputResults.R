#' Save the results
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param summarize if T a result averaged over all the replications is saved, if F each replication's result is saved
#'@param directory the directory to which the output will be saved (Enclose the name in double quotation!) (default: the current directory)
#'@param saveDataFileName the file name to save the simulated data with double-quotation, like "result1_1". (default: "BSLoutput")
#'
#'@return The simulation results (The output data was saved as BSLoutput.RData. After you load the data in R, you can find the data named as BSLoutput.)
#'
#'@export
outputResults <- function(sEnv=simEnv, summarize=T, directory=".", saveDataFileName="BSLoutput"){
  if(summarize){
    getMean <- function(data){
      tapply(data$gValue, data$genoRec$basePopID, mean)
    }
    getVar <- function(data){
      tapply(data$gValue, data$genoRec$basePopID, var)
    }
    muSim <- sapply(sEnv$sims, getMean)
    varSim <- sapply(sEnv$sims, getVar)
    BSLoutput <- cbind(muSim, varSim)
    colnames(BSLoutput) <- c(paste("mu", 1:sEnv$nSim, sep=""), paste("var", 1:sEnv$nSim, sep=""))
  }else{
    BSLoutput <- sEnv$sims
  }
  saveRDS(BSLoutput, file=paste(directory, "/", saveDataFileName, ".rds", sep=""))
}
