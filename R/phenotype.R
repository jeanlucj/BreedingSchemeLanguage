#'Evaluate the phenotypic value
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param errorVar error variance
#'@param popID population ID to be evaluated (default: the latest population)
#'@param locations integer vector of the locations where phenotyping occurs (e.g., c(1, 3) to phenotype at locations 1 and 3. Default: 1, phenotype at the first location)
#'@param years integer vector of the years when phenotyping occurs (e.g., 1:2 to phenotype during the first two years of the breeding scheme. Default: the year following any previous phenotyping. NOTE: thus, calling phenotype automatically increments the number of years the scheme takes)
#'
#'@return phenotypic values and the all information created before (list)
#'@export
# Locations and years get added when you phenotype in them for the first time
phenotype <- function(sEnv=simEnv, errorVar=1, popID=NULL, locations=1, years=ncol(sEnv$sims[[1]]$yearEffects)+1){
  parent.env(sEnv) <- environment()
  phenotype.func <- function(data, errorVar, popID, locations, years){
    # Who to phenotype
    if(is.null(popID)){
      popID <- max(data$genoRec$popID)
    }
    tf <- data$genoRec$popID %in% popID
    nPhen <- sum(tf)
    nLoc <- length(locations)
    nYr <- length(years)
    nRep <- nYr
    if (data$varParms$randLoc){
      nRep <- nRep * nLoc
      pValue <- calcPhenotypicValue(gv=data$gValue[tf,,drop=F], nRep, errorVar=errorVar)
    } else{
      if (!all(locations %in% 1:ncol(data$gValue))){
        stop("Phenotyping at unknown locations")
      }
      pValue <- calcPhenotypicValue(gv=data$gValue[tf, locations, drop=F], nRep, errorVar=errorVar)
    }
    # Year and location effects to add in
    nInd <- max(data$genoRec$GID)
    nAdd <- max(years) - ncol(data$yearEffects) # One col per year
    if (nAdd > 0 | (data$varParms$randLoc & max(locations) > ncol(data$locEffects))){
      M <- data$geno[1:nInd * 2 - 1, data$mapData$effectivePos] + data$geno[1:nInd * 2, data$mapData$effectivePos]
      nEffLoc <- length(data$mapData$effectivePos)
    }
    if (nAdd > 0){
      # Create GxY effects
      vp <- data$varParms$gByYearVar * data$varParms$fracGxEAdd
      gByYqtl <- matrix(rbinom(nEffLoc * nAdd, 1, 0.5), nEffLoc) * 2 - 1
      data$gByYqtl <- cbind(data$gByYqtl, gByYqtl)
      toAdd <- M %*% gByYqtl
      sdFound <- 1 / apply(toAdd[1:data$nFounders, , drop=F], 2, sd) * sqrt(vp)
      toAdd <- sapply(1:length(sdFound), function(i) toAdd[,i] * sdFound[i])
      data$yearScale <- c(data$yearScale, sdFound)
      data$yearEffects <- cbind(data$yearEffects, toAdd)
      vp <- data$varParms$gByYearVar * (1 - data$varParms$fracGxEAdd)
      toAdd <- matrix(rnorm(nInd * nAdd, sd=sqrt(vp)), nInd)
      data$yearEffectsI <- cbind(data$yearEffectsI, toAdd)
    }
    nAdd <- max(locations) - ncol(data$locEffects)
    if (data$varParms$randLoc & nAdd > 0){
      vp <- data$varParms$gByLocVar * data$varParms$fracGxEAdd
      gByLqtl <- matrix(rbinom(nEffLoc * nAdd, 1, 0.5), nEffLoc) * 2 - 1
      data$gByLqtl <- cbind(data$gByLqtl, gByLqtl)
      toAdd <- M %*% gByLqtl
      sdFound <- 1 / apply(toAdd[1:data$nFounders, , drop=F], 2, sd) * sqrt(vp)
      toAdd <- sapply(1:length(sdFound), function(i) toAdd[,i] * sdFound[i])
      data$locScale <- c(data$locScale, sdFound)
      data$locEffects <- cbind(data$locEffects, toAdd)
      vp <- data$varParms$gByLocVar * (1 - data$varParms$fracGxEAdd)
      toAdd <- matrix(rnorm(nInd * nAdd, sd=sqrt(vp)), nInd)
      data$locEffectsI <- cbind(data$locEffectsI, toAdd)
    }
    if (data$varParms$randLoc){
      pValue <- pValue + c(data$locEffects[tf, locations] + data$locEffectsI[tf, locations])
    }
    ye <- NULL
    for (i in years) ye <- c(ye, rep(data$yearEffects[tf,i], nLoc))
    pValue <- pValue + ye
    ye <- NULL
    for (i in years) ye <- c(ye, rep(data$yearEffectsI[tf,i], nLoc))
    pValue <- pValue + ye
    loc <- rep(locations, each=nPhen)
    yr <- rep(years, each=nPhen*nLoc)
    
    toAdd <- data.frame(phenoGID=data$genoRec$GID[tf], loc=loc, year=yr, error=errorVar, pValue=pValue)
    data$phenoRec <- rbind(data$phenoRec, toAdd)

    data$selCriterion <- list(popID=popID, criterion="pheno")
    # Take care of costs
    if (exists("totalCost", data)){
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
      sims <- sfLapply(sims, phenotype.func, errorVar=errorVar, popID=popID, locations=locations, years=years)
      sfStop()
    }else{
      sims <- lapply(sims, phenotype.func, errorVar=errorVar, popID=popID, locations=locations, years=years)
    }
  })
}
