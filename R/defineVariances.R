#' Define genetic, interaction, and error variances
#' 
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param gVariance genetic variance in the initial population Default is 1
#'@param locCorrelations matrix: genetic correlation in performance between locations default: NULL will cause locCorrelations to be gVariance / (gVariance + gByLocVar). If given, the genetic co-variance across locations is gVariance * locCorrelations. Default is NULL
#'@param gByLocVar scalar: the genotype by location variance Default is 1, BUT if locCorrelations given, this parameter is not used)
#'@param gByYearVar scalar: the genotype by year variance Default is 1
#'@param fracGxEAdd scalar: for GxL and GxY what fraction of the effect is additive versus non-additive Default is 0.8
#'@param plotTypeErrVars named vector: names are the plot types and contents are the error variances associated with them (default: Standard=1)
#'
#'@return Breeding scheme simulation object supplemented with variance parameters
#'
#'@export
defineVariances <- function(sEnv=simEnv, gVariance=1, locCorrelations=NULL, gByLocVar=1, gByYearVar=1, fracGxEAdd=0.8, plotTypeErrVars=c(Standard=1)){
  parent.env(sEnv) <- environment()
  variances.func <- function(bsl, gVariance, locCorrelations, gByLocVar, gByYearVar, fracGxEAdd, plotTypeErrVars){
    randLoc <- is.null(locCorrelations)
    if (randLoc){ # compound symmetric GxE here, with only g defined explicitly
      locCov <- matrix(gVariance)
    } else{ 
      locCov <- gVariance * locCorrelations
    }
    newEff <- nrow(locCov) - 1
    if (newEff > 0){
      bsl$mapData$effects <- cbind(bsl$mapData$effects, sapply(1:newEff, function(d) sample(bsl$mapData$effects[,1])))
    }
    
    bsl$varParms <- list(gVariance=gVariance, gByLocVar=gByLocVar, gByYearVar=gByYearVar, fracGxEAdd=fracGxEAdd, randLoc=randLoc, locCov=locCov, plotTypeErrVars=plotTypeErrVars)
    return(bsl)
  }
  
  with(sEnv, {
    # This is too fast to want to parallelize
    sims <- lapply(sims, variances.func, gVariance=gVariance, locCorrelations=locCorrelations, gByLocVar=gByLocVar, gByYearVar=gByYearVar, fracGxEAdd=fracGxEAdd, plotTypeErrVars=plotTypeErrVars)
  })
}
