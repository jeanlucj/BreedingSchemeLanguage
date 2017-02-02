#'Evaluate the phenotypic value
#'
#'@param errorVar error variance
#'@param popID population ID to be evaluated (default: the latest population)
#'
#'@return phenotypic values and the all information created before (list)
#'@export
phenotype <- function(simEnv, errorVar = 1, popID = NULL){
  parent.env(simEnv) <- environment()
  phenotype.func <- function(data, errorVar, popID){
    breedingData <- data$breedingData

    if(is.null(popID)){
      popID <- max(breedingData$popID)
    }
    tf <- breedingData$popID %in% popID
    
    pValue <- calcPhenotypicValue(gv = breedingData$gValue[tf], errorVar = errorVar)
    breedingData$pValue <- c(breedingData$pValue, pValue)
    breedingData$error <- c(breedingData$error, rep(errorVar, length(pValue)))
    breedingData$phenoGID <- c(breedingData$phenoGID, breedingData$GID[tf])
    
    data$breedingData <- breedingData
    data$selCriterion <- list(popID = popID, criterion = "pheno")
    return(data)
  }
  with(simEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, phenotype.func, errorVar = errorVar, popID = popID)
      sfStop()
    }else{
      sims <- lapply(sims, phenotype.func, errorVar = errorVar, popID = popID)
    }
  })
}
