#'Genotype markers
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param popID population ID to be genotyped (default: all populations)
#'
#'@return modifies the list sims in environment sEnv by indicating that genotypes of individuals in popID are available for breeding tasks
#'
#'@export
genotype <- function(sEnv=NULL, popID=NULL){
  genotype.func <- function(data, popID){
    nHasGeno <- sum(data$genoRec$hasGeno)
    if (is.null(popID)){
      data$genoRec$hasGeno <- TRUE
    } else{
      tf <- data$genoRec$popID %in% popID
      data$genoRec$hasGeno <- data$genoRec$hasGeno | tf
    }
    
    data$selCriterion$sharing <- "markers"
    
    if (exists("totalCost", data)) data$totalCost <- data$totalCost + (sum(data$genoRec$hasGeno) - nHasGeno) * data$costs$genoCost
    return(data)
  }
  
  if(is.null(sEnv)){
    if(exists("simEnv", .GlobalEnv)){
      sEnv <- get("simEnv", .GlobalEnv)
    } else{
      stop("No simulation environment was passed")
    }
  } 
  parent.env(sEnv) <- environment()
  with(sEnv, {
    # This is too fast to parallelize
    sims <- lapply(sims, genotype.func, popID=popID)
  })
}
