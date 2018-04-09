#' Save the results
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param summarize if T a result averaged over all the replications is saved, if F each replication's result is saved
#'@param saveDataFileName string of the file name to save the simulated data, like "result1_1". A path can be specified, like "simDirectory/result1_1" (in which case "simDirectory" must exist). Default: "BSLoutput".
#'
#'@return The output data is saved as a ".rds" file, with the file name dependent on the saveDataFileName parameter. To examine the data, use the readRDS function.
#'
#'@export
outputResults <- function(sEnv=NULL, summarize=T, saveDataFileName="BSLoutput"){
  if(is.null(sEnv)){
    if(exists("simEnv", .GlobalEnv)){
      sEnv <- get("simEnv", .GlobalEnv)
    } else{
      stop("No simulation environment was passed")
    }
  } 
  if(summarize){
    getMean <- function(data){
      tapply(data$gValue, data$genoRec$basePopID, mean)
    }
    getVar <- function(data){
      tapply(data$gValue, data$genoRec$basePopID, stats::var)
    }
    muSim <- sapply(sEnv$sims, getMean)
    varSim <- sapply(sEnv$sims, getVar)
    BSLoutput <- cbind(muSim, varSim)
    colnames(BSLoutput) <- c(paste("mu", 1:sEnv$nSim, sep=""), paste("var", 1:sEnv$nSim, sep=""))
  }else{
    BSLoutput <- sEnv$sims
  }
  saveRDS(BSLoutput, file=paste(saveDataFileName, ".rds", sep=""))
}
