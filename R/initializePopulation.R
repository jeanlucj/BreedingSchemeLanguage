#'Create a founder population
#'
#'@param simEnv an environment that BSL statements operate on
#'@param nInd population size
#'@param gVariance genetic variance in the initial population
#'
#'@return initial population informationand the all information created before (list)
#'
#'@export
initializePopulation <- function(simEnv, nInd=100, gVariance=1){
  parent.env(simEnv) <- environment()
  initializePopulation.func <- function(data, nInd, gVariance){
    mapData <- data$mapData
    founderHaps <- data$founderHaps
    seed <- round(runif(1, 0, 1e9))
    doubleGametes <- function(gametes){
      genotypes <- matrix(NA, 2 * nrow(gametes), ncol(gametes))
      genotypes[1:nrow(gametes) * 2, ] <- gametes
      genotypes[1:nrow(gametes) * 2 - 1, ] <- gametes
    }
    geno <- doubleGametes(founderHaps)
    geno <- randomMate(popSize=nrow(geno) / 2, geno=geno, pos=mapData$map$Pos)$progenies
    geno <- randomMate(popSize=nInd, geno=geno, pos=mapData$map$Pos)$progenies
    geno <- geno * 2 - 1
    gValue <- calcGenotypicValue(geno=geno, mapData=mapData)
    coef <- sqrt(gVariance / var(gValue))
    mapData$effects <- as.matrix(mapData$effects * coef, ncol=1)
    gValue <- coef*gValue
    GID <- 1:nInd
    popID <- rep(0, nInd)
    hasGeno <- rep(FALSE, nInd)
    breedingData <- list(geno=geno, GID=GID, popID=popID, popIDsel=popID, hasGeno=hasGeno, gValue=gValue)
    return(list(mapData=mapData, breedingData=breedingData))
  }
  with(simEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, initializePopulation.func, nInd=nInd, gVariance=gVariance)
      sfStop()
    }else{
      sims <- lapply(sims, initializePopulation.func, nInd=nInd, gVariance=gVariance)
    }
  })
}
