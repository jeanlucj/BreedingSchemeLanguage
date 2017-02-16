#'Evaluate the phenotypic value
#'
#'@param simEnv an environment that BSL statements operate on
#'@param errorVar error variance
#'@param popID population ID to be evaluated (default: the latest population)
#'
#'@return phenotypic values and the all information created before (list)
#'@export
phenotype <- function(simEnv, errorVar=1, popID=NULL){
  parent.env(simEnv) <- environment()
  phenotype.func <- function(data, errorVar, popID){

    if(is.null(popID)){
      popID <- max(data$genoRec$popID)
    }
    tf <- data$genoRec$popID %in% popID
    
    pValue <- calcPhenotypicValue(gv=data$genoRec$gValue[tf], errorVar=errorVar)
    
    toAdd <- data.frame(phenoGID=data$genoRec$GID[tf], pValue=pValue, error=errorVar)
    data$phenoRec <- rbind(data$phenoRec, toAdd)

    data$selCriterion <- list(popID=popID, criterion="pheno")
    return(data)
  }
  with(simEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, phenotype.func, errorVar=errorVar, popID=popID)
      sfStop()
    }else{
      sims <- lapply(sims, phenotype.func, errorVar=errorVar, popID=popID)
    }
  })
}
