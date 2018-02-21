#'Evaluate the phenotypic value
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param plotType links to previously defined error variance and cost (default: "Standard")
#'@param nRep scalar: the number of replications per trial (i.e., year x loc combination) (default: 1)
#'@param popID population ID to be evaluated (default: the latest population)
#'@param locations integer vector of the locations where phenotyping occurs (e.g., c(1, 3) to phenotype at locations 1 and 3. Default: 1, phenotype at the first location)
#'@param years integer vector of the years when phenotyping occurs (e.g., 1:2 to phenotype during the first two years of the breeding scheme. Default: the last year among previous phenotyping. NOTE: thus, to phenotype in a new [the next] year, specify the next year number [e.g., if past phenotyping was in years 1 & 2, specify 3]).
#'
#'@return modifies the list sims in environment sEnv by generating phenotypes for the specified popID, locations, and years.
#'
#'@export
# Locations and years get added when you phenotype in them for the first time
phenotype <- function(sEnv=NULL, plotType="Standard", nRep=1, popID=NULL, locations=1, years=NULL){
  phenotype.func <- function(bsl, plotType, nRep, popID, locations, years){
    errorVar <- bsl$varParms$plotTypeErrVars[plotType] / nRep
    names(errorVar) <- NULL
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
    nTrial <- nYr
    if (bsl$varParms$randLoc){
      nTrial <- nTrial * nLoc
      pValue <- calcPhenotypicValue(gv=bsl$gValue[tf,,drop=F], nRep=nTrial, errorVar=errorVar)
    } else{
      if (!all(locations %in% 1:ncol(bsl$gValue))){
        stop("Phenotyping at unknown locations")
      }
      pValue <- calcPhenotypicValue(gv=bsl$gValue[tf, locations, drop=F], nRep=nTrial, errorVar=errorVar)
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
      gByYqtl <- matrix(stats::rbinom(nEffLoc * nAdd, 1, 0.5), nEffLoc) * 2 - 1
      bsl$gByYqtl <- cbind(bsl$gByYqtl, gByYqtl)
      toAdd <- M %*% gByYqtl
      sdFound <- 1 / apply(toAdd[1:bsl$nFounders, , drop=F], 2, stats::sd) * sqrt(vp)
      toAdd <- sapply(1:length(sdFound), function(i) toAdd[,i] * sdFound[i])
      bsl$yearScale <- c(bsl$yearScale, sdFound)
      bsl$yearEffects <- cbind(bsl$yearEffects, toAdd)
      vp <- bsl$varParms$gByYearVar * (1 - bsl$varParms$fracGxEAdd)
      toAdd <- matrix(stats::rnorm(nInd * nAdd, sd=sqrt(vp)), nInd)
      bsl$yearEffectsI <- cbind(bsl$yearEffectsI, toAdd)
      if (exists("totalCost", bsl)){
        bsl$totalCost <- bsl$totalCost + nAdd * bsl$costs$yearCost
      }
    }
    nAdd <- max(locations) - ncol(bsl$locEffects)
    if (bsl$varParms$randLoc & nAdd > 0){
      vp <- bsl$varParms$gByLocVar * bsl$varParms$fracGxEAdd
      gByLqtl <- matrix(stats::rbinom(nEffLoc * nAdd, 1, 0.5), nEffLoc) * 2 - 1
      bsl$gByLqtl <- cbind(bsl$gByLqtl, gByLqtl)
      toAdd <- M %*% gByLqtl
      sdFound <- 1 / apply(toAdd[1:bsl$nFounders, , drop=F], 2, stats::sd) * sqrt(vp)
      toAdd <- sapply(1:length(sdFound), function(i) toAdd[,i] * sdFound[i])
      bsl$locScale <- c(bsl$locScale, sdFound)
      bsl$locEffects <- cbind(bsl$locEffects, toAdd)
      vp <- bsl$varParms$gByLocVar * (1 - bsl$varParms$fracGxEAdd)
      toAdd <- matrix(stats::rnorm(nInd * nAdd, sd=sqrt(vp)), nInd)
      bsl$locEffectsI <- cbind(bsl$locEffectsI, toAdd)
      if (exists("totalCost", bsl)){
        bsl$totalCost <- bsl$totalCost + nAdd * bsl$costs$locCost
      }
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

    bsl$selCriterion <- list(popID=popID, criterion="pheno", sharing="none")
    # Take care of costs
    if (exists("totalCost", bsl)){
      perPlotCost <- bsl$costs$phenoCost[plotType]
      bsl$totalCost <- bsl$totalCost + nPhen * perPlotCost * nLoc * nYr * nRep
    }
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
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, phenotype.func, plotType=plotType, nRep=nRep, popID=popID, locations=locations, years=years)
      sfStop()
    }else{
      sims <- lapply(sims, phenotype.func, plotType=plotType, nRep=nRep, popID=popID, locations=locations, years=years)
    }
  })
}
