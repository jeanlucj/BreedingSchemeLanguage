#' Save the results
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param summarize If TRUE a matrix with breeding cycles in rows and replications in columns, with the first set of columns being cycle means and the second set of columns cycle variances. If FALSE a list as long as the number of replications, with each list element containing a list of all simulation objects.
#'@param saveDataFileName NULL or string of the file name to save the simulated data, like "result1_1". A path can be specified, like "simDirectory/result1_1" (in which case "simDirectory" must exist). Default: NULL.
#'
#'@return If saveDataFileName is NULL, data is return as an object, else data is saved as a ".rds" file. To examine the data, use the readRDS function.
#'
#'@export
outputResults <- function(sEnv=NULL, summarize=T, saveDataFileName=NULL){
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
    muSim <- matrix(sapply(sEnv$sims, getMean), ncol=sEnv$nSim)
    varSim <- matrix(sapply(sEnv$sims, getVar), ncol=sEnv$nSim)
    BSLoutput <- cbind(muSim, varSim)
    colnames(BSLoutput) <- c(paste("mu", 1:sEnv$nSim, sep=""), paste("var", 1:sEnv$nSim, sep=""))
  }else{
    BSLoutput <- sEnv$sims
  }
  if (is.null(saveDataFileName)) return(BSLoutput) else saveRDS(BSLoutput, file=paste(saveDataFileName, ".rds", sep=""))
}
