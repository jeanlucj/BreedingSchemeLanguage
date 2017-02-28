#'Evaluate the phenotypic value
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param errorVar error variance
#'@param popID population ID to be evaluated (default: the latest population)
#'@param locations integer vector of the locations where phenotyping occurs (e.g., c(1, 3) to phenotype at locations 1 and 3. Default: 1, phenotype at the first location)
#'@param years integer vector of the years when phenotyping occurs (e.g., 1:2 to phenotype during the first two years of the breeding scheme. Default: the last year among previous phenotyping. NOTE: thus, to phenotype in a new [the next] year, specify the next year number [e.g., if past phenotyping was in years 1 & 2, specify 3]).
#'
#'@return phenotypic values and the all information created before (list)
#'@export
# Locations and years get added when you phenotype in them for the first time
phenotype <- function(sEnv=simEnv, errorVar=1, popID=NULL, locations=1, years=NULL){
  parent.env(sEnv) <- environment()
  phenotype.func <- function(bsl, errorVar, popID, locations, years){
    # When to phenotype
    if (is.null(years)) years=max(ncol(bsl$yearEffects), 1)
    # Who to phenotype
    if(is.null(popID)){
      popID <- max(bsl$genoRec$popID)
    }
    tf <- bsl$genoRec$popID %in% popID
    nPhen <- sum(tf)
    nLoc <- length(locations)
    nYr <- length(years)
    nRep <- nYr
    if (bsl$varParms$randLoc){
      nRep <- nRep * nLoc
      pValue <- calcPhenotypicValue(gv=bsl$gValue[tf,,drop=F], nRep, errorVar=errorVar)
    } else{
      if (!all(locations %in% 1:ncol(bsl$gValue))){
        stop("Phenotyping at unknown locations")
      }
      pValue <- calcPhenotypicValue(gv=bsl$gValue[tf, locations, drop=F], nRep, errorVar=errorVar)
    }
    # Year and location effects to add in
    nInd <- max(bsl$genoRec$GID)
    nAdd <- max(years) - ncol(bsl$yearEffects) # One col per year
    if (nAdd > 0 | (bsl$varParms$randLoc & max(locations) > ncol(bsl$locEffects))){
      M <- bsl$geno[1:nInd * 2 - 1, bsl$mapData$effectivePos] + bsl$geno[1:nInd * 2, bsl$mapData$effectivePos]
      nEffLoc <- length(bsl$mapData$effectivePos)
    }
    if (nAdd > 0){
      # Create GxY effects
      vp <- bsl$varParms$gByYearVar * bsl$varParms$fracGxEAdd
      gByYqtl <- matrix(rbinom(nEffLoc * nAdd, 1, 0.5), nEffLoc) * 2 - 1
      bsl$gByYqtl <- cbind(bsl$gByYqtl, gByYqtl)
      toAdd <- M %*% gByYqtl
      sdFound <- 1 / apply(toAdd[1:bsl$nFounders, , drop=F], 2, sd) * sqrt(vp)
      toAdd <- sapply(1:length(sdFound), function(i) toAdd[,i] * sdFound[i])
      bsl$yearScale <- c(bsl$yearScale, sdFound)
      bsl$yearEffects <- cbind(bsl$yearEffects, toAdd)
      vp <- bsl$varParms$gByYearVar * (1 - bsl$varParms$fracGxEAdd)
      toAdd <- matrix(rnorm(nInd * nAdd, sd=sqrt(vp)), nInd)
      bsl$yearEffectsI <- cbind(bsl$yearEffectsI, toAdd)
    }
    nAdd <- max(locations) - ncol(bsl$locEffects)
    if (bsl$varParms$randLoc & nAdd > 0){
      vp <- bsl$varParms$gByLocVar * bsl$varParms$fracGxEAdd
      gByLqtl <- matrix(rbinom(nEffLoc * nAdd, 1, 0.5), nEffLoc) * 2 - 1
      bsl$gByLqtl <- cbind(bsl$gByLqtl, gByLqtl)
      toAdd <- M %*% gByLqtl
      sdFound <- 1 / apply(toAdd[1:bsl$nFounders, , drop=F], 2, sd) * sqrt(vp)
      toAdd <- sapply(1:length(sdFound), function(i) toAdd[,i] * sdFound[i])
      bsl$locScale <- c(bsl$locScale, sdFound)
      bsl$locEffects <- cbind(bsl$locEffects, toAdd)
      vp <- bsl$varParms$gByLocVar * (1 - bsl$varParms$fracGxEAdd)
      toAdd <- matrix(rnorm(nInd * nAdd, sd=sqrt(vp)), nInd)
      bsl$locEffectsI <- cbind(bsl$locEffectsI, toAdd)
    }
    if (bsl$varParms$randLoc){
      pValue <- pValue + c(bsl$locEffects[tf, locations] + bsl$locEffectsI[tf, locations])
    }
    ye <- NULL
    for (i in years) ye <- c(ye, rep(bsl$yearEffects[tf,i], nLoc))
    pValue <- pValue + ye
    ye <- NULL
    for (i in years) ye <- c(ye, rep(bsl$yearEffectsI[tf,i], nLoc))
    pValue <- pValue + ye
    loc <- rep(locations, each=nPhen)
    yr <- rep(years, each=nPhen*nLoc)
    
    toAdd <- data.frame(phenoGID=bsl$genoRec$GID[tf], loc=loc, year=yr, error=errorVar, pValue=pValue)
    bsl$phenoRec <- rbind(bsl$phenoRec, toAdd)

    bsl$selCriterion <- list(popID=popID, criterion="pheno")
    # Take care of costs
    if (exists("totalCost", bsl)){
      perPlotCost <- abs(bsl$costs$phenoCost$error - errorVar)
      perPlotCost <- which(perPlotCost == min(perPlotCost))
      perPlotCost <- bsl$costs$phenoCost$cost[perPlotCost]
      bsl$totalCost <- bsl$totalCost + nPhen * perPlotCost * nLoc * nYr
    }
    return(bsl)
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
