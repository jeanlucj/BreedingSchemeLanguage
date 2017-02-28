#' Define the relationships among locations
#' NOTE: this function needs to be called after defineSpecies
#' 
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param gVariance genetic variance in the initial population Default is 1
#'@param gByLocVar scalar: the genotype by location variance Default is 1, BUT if locCorrelations given, this parameter is not used)
#'@param gByYearVar scalar: the genotype by year variance Default is 1
#'@param fracGxEAdd scalar: for GxL and GxY what fraction of the effect is additive versus non-additive Default is 0.8
#'@param locCorrelations matrix: genetic correlation in performance between locations default: NULL will cause locCorrelations to be gVariance / (gVariance + gByLocVar). If given, the genetic co-variance across locations is gVariance * locCorrelations. Default is NULL
#'
#'@return Species information and input values for the simulation (list)
#'
#'@export
defineVariances <- function(sEnv=simEnv, gVariance=1, locCorrelations=NULL, gByLocVar=1, gByYearVar=1, fracGxEAdd=0.8){
  parent.env(sEnv) <- environment()
  variances.func <- function(bsl, gVariance, locCorrelations, gByLocVar, gByYearVar, fracGxEAdd){
    randLoc <- is.null(locCorrelations)
    if (randLoc){ # compound symmetric GxE here, with only g defined explicitly
      locCov <- matrix(gVariance)
    } else{ 
      locCov <- gVariance * locCorrelations
    }
    nQTL <- max(bsl$mapData$effectID)
    bsl$mapData$effects <- matrix(rnorm(nQTL * nrow(locCov)), nQTL) %*% chol(locCov)
    bsl$varParms <- list(gVariance=gVariance, gByLocVar=gByLocVar, gByYearVar=gByYearVar, fracGxEAdd=fracGxEAdd, randLoc=randLoc, locCov=locCov)
    return(bsl)
  }
  with(sEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, variances.func, gVariance=gVariance, locCorrelations=locCorrelations, gByLocVar=gByLocVar, gByYearVar=gByYearVar, fracGxEAdd=fracGxEAdd)
      sfStop()
    }else{
      sims <- lapply(sims, variances.func, gVariance=gVariance, locCorrelations=locCorrelations, gByLocVar=gByLocVar, gByYearVar=gByYearVar, fracGxEAdd=fracGxEAdd)
    }
  })
}
