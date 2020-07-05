#'Genotype markers
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param popID population ID to be genotyped (default: all populations)
#'
#'@return modifies the list sims in environment sEnv by indicating that genotypes of individuals in popID are available for breeding tasks
#'
#'@export
genotype <- function(sEnv=NULL, popID=NULL){
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
      hasGenoBefore <- sum(budgetRec$hasGeno)
      if (is.null(popID)){
        budgetRec$hasGeno <- TRUE
      } else{
        budgetRec$hasGeno <- budgetRec$hasGeno | (budgetRec$popID %in% popID)
      }
      totalCost <- totalCost + (sum(budgetRec$hasGeno) - hasGenoBefore) * costs$genoCost
      rm(hasGenoBefore)
    }
    # This is too fast to parallelize
    if (!exists(onlyCost)) onlyCost <- FALSE
    if (!onlyCost){
      sims <- lapply(sims, genotype.func, popID=popID)
    }
  })
}
