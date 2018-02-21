#'Self-fertilize
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param nProgeny the number of progeny
#'@param popID population ID to be self-fertilized (default: the latest population)
#'
#'@return modifies the list sims in environment sEnv by creating a selfed progeny population as specified, with an incremented population number
#'
#'@export
selfFertilize <- function(sEnv=NULL, nProgeny=100, popID=NULL){
  selfFertilize.func <- function(data, nProgeny, popID){
    locPos <- data$mapData$map$Pos
    if(is.null(popID)){
      popID <- max(data$genoRec$popID)
    }
    tf <- data$genoRec$popID %in% popID
    GIDpar <- data$genoRec$GID[tf]
    nPar <- length(GIDpar)
    geno <- data$geno[rep(GIDpar*2, each=2) + rep(-1:0, nPar),]
    geno <- makeSelfs(popSize=nProgeny, geno=geno, pos=locPos)
    pedigree <- cbind(matrix(GIDpar[geno$pedigree], nProgeny), 1)
    geno <- geno$progenies
    data <- addProgenyData(data, geno, pedigree)
    if (exists("totalCost", data)) data$totalCost <- data$totalCost + nProgeny * data$costs$selfCost
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
    if(nCore > 1){
      snowfall::sfInit(parallel=T, cpus=nCore)
      sims <- snowfall::sfLapply(sims, selfFertilize.func, nProgeny=nProgeny, popID=popID)
      snowfall::sfStop()
    }else{
      sims <- lapply(sims, selfFertilize.func, nProgeny=nProgeny, popID=popID)
    }
  })
}
