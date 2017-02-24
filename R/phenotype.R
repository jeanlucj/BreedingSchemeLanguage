#'Evaluate the phenotypic value
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param errorVar error variance
#'@param popID population ID to be evaluated (default: the latest population)
#'@param locations integer vector of the locations where phenotyping occurs (e.g., c(1, 3) to phenotype at locations 1 and 3. Default: 1, phenotype at the first location)
#'@param years integer vector of the years when phenotypin occurs (e.g., 1:2 to phenotype during the first two years of the breeding scheme. Default: 1, phenotype the first year)
#'
#'@return phenotypic values and the all information created before (list)
#'@export
# Locations and years get added when you phenotype in them for the first time
phenotype <- function(sEnv=simEnv, errorVar=1, popID=NULL, locations=1, years=1){
  parent.env(sEnv) <- environment()
  phenotype.func <- function(data, errorVar, popID, locations, years){
    # Who to phenotype
    if(is.null(popID)){
      popID <- max(data$genoRec$popID)
    }
    tf <- data$genoRec$popID %in% popID
    pValue <- calcPhenotypicValue(gv=data$gValue[tf,], errorVar=errorVar)
    nPhen <- nrow(pValue)
    
    # Year and location effects to add in
    nInd <- max(data$genoRec$GID)
    nAdd <- max(years) - ncol(data$yearEffects)
    if (nAdd > 0){
      vp <- data$varParms$gByYearVar * data$varParms$fracGxEAdd
      toAdd <- rmvnorm(nAdd, mu=rep(0, nInd), Sigma=data$qtlRelMat) * sqrt(vp)
      data$yearEffects <- cbind(data$yearEffects, toAdd)
      vp <- data$varParms$gByYearVar * (1 - data$varParms$fracGxEAdd)
      toAdd <- matrix(rnorm(nInd * nAdd, sd=sqrt(vp)), nInd)
      data$yearEffectsI <- cbind(data$yearEffectsI, toAdd)
    }
    if (data$randLoc){
      nAdd <- max(locations) - ncol(data$locEffects)
      if (nAdd > 0){
        vp <- data$varParms$gByLocVar * data$varParms$fracGxEAdd
        toAdd <- rmvnorm(nAdd, mu=rep(0, nInd), Sigma=data$qtlRelMat) * sqrt(vp)
        data$locEffects <- cbind(data$locEffects, toAdd)
        vp <- data$varParms$gByLocVar * (1 - data$varParms$fracGxEAdd)
        toAdd <- matrix(rnorm(nInd * nAdd, sd=sqrt(vp)), nInd)
        data$locEffectsI <- cbind(data$locEffectsI, toAdd)
      }
    }
    nLoc <- length(locations)
    nYr <- length(years)
    if (data$randLoc){
      pValue <- rep(pValue, nLoc * nYr)
      pValue <- pValue + data$locEffects[tf, locations] + data$locEffectsI[tf, locations]
    } else{
      pValue <- rep(pValue[tf, locations], nYr)
    }
    ye <- NULL; for (i in years) ye <- c(ye, rep(data$yearEffects[,i], nLoc))
    pValue <- pValue + ye
    ye <- NULL; for (i in years) ye <- c(ye, rep(data$yearEffectsI[,i], nLoc))
    pValue <- pValue + ye
    locFact <- as.factor(rep(locations, each=nPhen))
    yrFact <- as.factor(rep(years, each=nPhen*nLoc))
    
    toAdd <- data.frame(phenoGID=as.factor(data$genoRec$GID[tf]), loc=locFact, year=yrFact, error=errorVar, pValue=pValue)
    data$phenoRec <- rbind(data$phenoRec, toAdd)

    data$selCriterion <- list(popID=popID, criterion="pheno")
    # Take care of costs
    if (exists("totalCost", data){
      perPlotCost <- abs(data$costs$phenoCost$error - errorVar)
      perPlotCost <- which(perPlotCost == min(perPlotCost))
      perPlotCost <- data$costs$phenoCost$cost[perPlotCost]
      data$totalCost <- data$totalCost + nPhen * perPlotCost * nLoc * nYr
    }
    return(data)
  }
  with(sEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, phenotype.func, errorVar=errorVar, popID=popID)
      sfStop()
    }else{
      sims <- lapply(sims, phenotype.func, errorVar=errorVar, popID=popID)
    }
  })
}
