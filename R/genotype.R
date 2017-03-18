#'Genotype markers
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param popID population ID to be genotyped (default: all populations)
#'@return marker genotype and the all information created before (list)
#'
#'@export
genotype <- function(sEnv=simEnv, popID=NULL){
  parent.env(sEnv) <- environment()
  genotype.func <- function(data, popID){
    nHasGeno <- sum(data$genoRec$hasGeno)
    if (is.null(popID)){
      data$genoRec$hasGeno <- TRUE
    } else{
      tf <- data$genoRec$popID %in% popID
      data$genoRec$hasGeno <- data$genoRec$hasGeno | tf
    }
    if (exists("totalCost", data)) data$totalCost <- data$totalCost + (sum(data$genoRec$hasGeno) - nHasGeno) * data$costs$genoCost
    return(data)
  }
  with(sEnv, {
    # This is too fast to want to parallelize
    sims <- lapply(sims, genotype.func, popID=popID)
  })
}
