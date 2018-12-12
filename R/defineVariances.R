#' Define genetic, interaction, and error variances
#' 
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param gVariance genetic variance in the initial population Default is 1
#'@param locCorrelations matrix: genetic correlation in performance between locations default: If given, the genetic co-variance across locations is gVariance * locCorrelations. If NULL, deviations with variance gByLocVar will be added to the main genotypic effect to determine the genotypic value in each location. The expected genetic correlation between locations is then gVariance / (gVariance + gByLocVar). Default is NULL
#'@param gByLocVar scalar: the genotype by location variance. Default is 1, BUT if locCorrelations given, this parameter is not used.
#'@param gByYearVar scalar: the genotype by year variance. Default is 1.
#'@param fracGxEAdd scalar: for GxL and GxY what fraction of the effect is additive versus non-additive. Default is 1. NOTE: the additive component of GxL and GxY is due to the QTL generating the main genotypic effect. Non-additive components of GxL and GxY are sampled at random, normally distributed, IID.
#'@param plotTypeErrVars named vector: names are the plot types and contents are the error variances associated with them (default: Standard=1)
#'
#'@return modifies the list sims in environment sEnv by adding parameters that determine genetic, location, year, and error variances
#'
#'@export
defineVariances <- function(sEnv=NULL, gVariance=1, locCorrelations=NULL, gByLocVar=1, gByYearVar=1, fracGxEAdd=1, plotTypeErrVars=c(Standard=1)){
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
  
  if(is.null(sEnv)){
    if(exists("simEnv", .GlobalEnv)){
      sEnv <- get("simEnv", .GlobalEnv)
    } else{
      stop("No simulation environment was passed")
    }
  } 
  parent.env(sEnv) <- environment()
  with(sEnv, {
    # This is too fast to parallelize
    sims <- lapply(sims, variances.func, gVariance=gVariance, locCorrelations=locCorrelations, gByLocVar=gByLocVar, gByYearVar=gByYearVar, fracGxEAdd=fracGxEAdd, plotTypeErrVars=plotTypeErrVars)
  })
}
