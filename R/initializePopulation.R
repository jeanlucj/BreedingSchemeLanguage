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
    seed <- round(runif(1, 0, 1e9))
    mapData <- data$mapData
    geno <- data$founderHaps * 2 - 1
    geno <- geno[sample(nrow(geno), nInd*2, replace=T),]
    geno <- randomMate(popSize=nInd, geno=geno, pos=mapData$map$Pos)
    pedigree <- -geno$pedigree # For founders, parents will be negative
    geno <- geno$progenies
    gValue <- calcGenotypicValue(geno=geno, mapData=mapData)
    coef <- sqrt(gVariance / var(gValue))
    mapData$effects <- as.matrix(mapData$effects * coef, ncol=1)
    gValue <- coef*gValue
    GID <- 1:nInd
    popID <- rep(0, nInd)
    hasGeno <- rep(FALSE, nInd)
    genoRec <- data.frame(GID=GID, pedigree=pedigree, popID=popID, basePopID=popID, hasGeno=hasGeno, gValue=gValue)
    return(list(mapData=mapData, geno=geno, genoRec=genoRec))
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
