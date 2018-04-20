#'Genotype markers
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param popID population ID to be genotyped (default: all populations)
#'@param parms an optional named list or vector. Objects with those names will be created with the corresponding values. A way to pass values that are not predetermined by the script.
#'
#'@return modifies the list sims in environment sEnv by indicating that genotypes of individuals in popID are available for breeding tasks
#'
#'@export
genotype <- function(sEnv=NULL, popID=NULL, parms=NULL){
  if(!is.null(parms)){
    for (n in 1:length(parms)){
      assign(names(parms)[n], parms[[n]])
    }
  }
  genotype.func <- function(bsl, popID){
    nHasGeno <- sum(bsl$genoRec$hasGeno)
    if (is.null(popID)){
      bsl$genoRec$hasGeno <- TRUE
    } else{
      tf <- bsl$genoRec$popID %in% popID
      bsl$genoRec$hasGeno <- bsl$genoRec$hasGeno | tf
    }
    
    bsl$selCriterion$sharing <- "markers"
    
    if (exists("totalCost", bsl)) bsl$totalCost <- bsl$totalCost + (sum(bsl$genoRec$hasGeno) - nHasGeno) * bsl$costs$genoCost
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
    # This is too fast to parallelize
    sims <- lapply(sims, genotype.func, popID=popID)
  })
}
