#'Self-fertilize
#'
#'@param simEnv an environment that BSL statements operate on
#'@param nProgeny the number of progeny
#'@param popID population ID to be self-fertilized (default: the latest population)
#'
#'@return sequence information of progenies and the all information created before (list)
#'
#'@export
selfFertilize <- function(simEnv, nProgeny=100, popID=NULL){
  parent.env(simEnv) <- environment()
  selfFertilize.func <- function(data, nProgeny, popID){
    locPos <- data$mapData$map$Pos
    breedingData <- data$breedingData
    if(is.null(popID)){
      popID <- max(breedingData$popID)
    }
    tf <- breedingData$popID %in% popID
    GID.now <- breedingData$GID[tf]
    geno.now <- breedingData$geno[sort(c(GID.now * 2 - 1, GID.now * 2)), ]
    geno.progeny <- makeSelfs(popSize=nProgeny, geno=geno.now, pos=locPos)$progenies
    gValue <- calcGenotypicValue(geno=geno.progeny, mapData=data$mapData)
    GID.progeny <- max(breedingData$GID) + 1:(nrow(geno.progeny) / 2)
    breedingData$GID <- c(breedingData$GID, GID.progeny)
    popID.progeny <- rep(max(breedingData$popID) + 1, nrow(geno.progeny) / 2)
    breedingData$popID <- c(breedingData$popID, popID.progeny)
    breedingData$popIDsel <- c(breedingData$popIDsel, popID.progeny)
    breedingData$geno <- rbind(breedingData$geno, geno.progeny)
    breedingData$hasGeno <- c(breedingData$hasGeno, rep(FALSE, nProgeny))
    breedingData$gValue <- c(breedingData$gValue, gValue)
    data$breedingData <- breedingData
    return(data)
  }
  with(simEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, selfFertilize.func, nProgeny=nProgeny, popID=popID)
      sfStop()
    }else{
      sims <- lapply(sims, selfFertilize.func, nProgeny=nProgeny, popID=popID)
    }
  })
}
