#'Genotype markers
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param popID population ID to be genotyped (default: all populations)
#'@param onlyCost logical. If true, don't do the breeding task, just calculate its cost.  Default: FALSE. 
#'
#'@return modifies the list sims in environment sEnv by indicating that genotypes of individuals in popID are available for breeding tasks
#'
#'@export
genotype <- function(sEnv=NULL, popID=NULL, onlyCost=FALSE){
  genotype.func <- function(bsl, popID){
    nHasGeno <- sum(bsl$genoRec$hasGeno)
    if (is.null(popID)){
      bsl$genoRec$hasGeno <- TRUE
    } else{
      tf <- bsl$genoRec$popID %in% popID
      bsl$genoRec$hasGeno <- bsl$genoRec$hasGeno | tf
    }
    bsl$selCriterion$sharing <- "markers"
    
    return(bsl)
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
    if (exists("totalCost")){
      genoRec <- sims[[1]]$genoRec
      nHasGeno <- sum(genoRec$hasGeno)
      if (is.null(popID)){
        genoRec$hasGeno <- TRUE
      } else{
        tf <- genoRec$popID %in% popID
        genoRec$hasGeno <- genoRec$hasGeno | tf
      }
      totalCost <- totalCost + (sum(genoRec$hasGeno) - nHasGeno) * costs$genoCost
    }
    # This is too fast to parallelize
    if (!onlyCost){
      sims <- lapply(sims, genotype.func, popID=popID)
    }
  })
}
