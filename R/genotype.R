#'Genotype markers
#'
#'@param simEnv an environment that BSL statements operate on
#'@param popID population ID to be genotyped (default: all populations)
#'@return marker genotype and the all information created before (list)
#'
#'@export
genotype <- function(simEnv, popID=NULL){
  parent.env(simEnv) <- environment()
  genotype.func <- function(data){
    if (is.null(popID)){
      data$breedingData$hasGeno <- rep(TRUE, length(data$breedingData$GID))
    } else{
      tf <- data$breedingData$popID %in% popID
      data$breedingData$hasGeno <- data$breedingData$hasGeno | tf
    }
    return(data)
  }
  with(simEnv, {
    # This is too fast to want to parallelize
    sims <- lapply(sims, genotype.func)
  })
}
