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
  selfFertilize.func <- function(bsl, nProgeny, popID){
    locPos <- bsl$mapData$map$Pos
    if(is.null(popID)){
      popID <- max(bsl$genoRec$popID)
    }
    tf <- bsl$genoRec$popID %in% popID
    GIDpar <- bsl$genoRec$GID[tf]
    nPar <- length(GIDpar)
    geno <- bsl$geno[rep(GIDpar*2, each=2) + rep(-1:0, nPar),]
    geno <- makeSelfs(popSize=nProgeny, geno=geno, pos=locPos)
    pedigree <- cbind(matrix(GIDpar[geno$pedigree], nProgeny), 1)
    geno <- geno$progenies
    bsl <- addProgenyData(bsl, geno, pedigree)
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
      budgetRec <- rbind(budgetRec, data.frame(GID=max(budgetRec$GID) + 1:nProgeny, popID=rep(max(budgetRec$popID) + 1, nProgeny), hasGeno=FALSE))
      totalCost <- totalCost + nProgeny * costs$selfCost 
    }
    
    if (!onlyCost){
      if(nCore > 1){
        snowfall::sfInit(parallel=T, cpus=nCore)
        sims <- snowfall::sfLapply(sims, selfFertilize.func, nProgeny=nProgeny, popID=popID)
        snowfall::sfStop()
      }else{
        sims <- lapply(sims, selfFertilize.func, nProgeny=nProgeny, popID=popID)
      }
    }
  })
}
