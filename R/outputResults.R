#' Save the results
#'
#'@param simEnv an environment that BSL statements operate on
#'@param summarize if T a result averaged over all the replications is saved, if F each replication's result is saved
#'@param directory the directory to which the output will be saved (Enclose the name in double quotation!) (default: the current directory)
#'@param saveDataFileName the file name to save the simulated data with double-quotation, like "result1_1". (default: "BSLoutput")
#'
#'@return The simulation results (The output data was saved as BSLoutput.RData. After you load the data in R, you can find the data named as BSLoutput.)
#'
#'@export
outputResults <- function(simEnv, summarize=T, directory=".", saveDataFileName="BSLoutput"){
  if(summarize){
    getMean <- function(data){
      tapply(data$genoRec$gValue, data$genoRec$basePopID, mean)
    }
    getVar <- function(data){
      tapply(data$genoRec$gValue, data$genoRec$basePopID, var)
    }
    muSim <- sapply(simEnv$sims, getMean)
    varSim <- sapply(simEnv$sims, getVar)
    BSLoutput <- cbind(muSim, varSim)
    # rownames(BSLoutput) <- popID
    colnames(BSLoutput) <- c(paste("mu", 1:simEnv$nSim, sep=""), paste("var", 1:simEnv$nSim, sep=""))
  }else{
    BSLoutput <- simEnv$sims
  }
  saveRDS(BSLoutput, file=paste(directory, "/", saveDataFileName, ".rds", sep=""))
}
