#' Define the relationships among locations
#' NOTE: if this function is not called, genotype by location and genotype by year effects
#' are not used in the simulation
#' 
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param gVariance genetic variance in the initial population
#'@param gByLocVar scalar: the genotype by location variance (default: 1 BUT if locCorrelations given, this parameter is not used)
#'@param gByYearVar scalar: the genotype by year variance (default: 1)
#'@param fracGxEAdd scalar: for GxL and GxY what fraction of the effect is additive versus non-additive
#'@param nLoc scalar: number of locations default: if locCorrelations given, nrow(locCorrelations), else 2
#'@param locCorrelations matrix: genetic correlation in performance between locations default: NULL will cause locCorrelations to be gVariance / (gVariance + gByLocVar). If given, the genetic variance-covariance across locations is gVariance * locCorrelations
#'
#'@return Species information and input values for the simulation (list)
#'
#'@export
defineVariances <- function(sEnv=simEnv, gVariance=1, locCorrelations=NULL, nLoc=2, gByLocVar=1, gByYearVar=1, fracGxEAdd=0.8){
  parent.env(sEnv) <- environment()
  loc.func <- function(data){
    randLoc <- is.null(locCorrelations)
    if (randLoc){ # compound symmetric GxE here, with only g defined explicitly
      locCov <- gVariance
    } else{ 
      nLoc <- nrow(locCorrelations)
      locCov <- gVariance * locCorrelations
    }
    nQTL <- max(data$mapData$effectID)
    data$mapData$effects <- matrix(rnorm(nQTL * nrow(locCov)), nQTL) %*% chol(locCov)
    data$varParms <- list(gVariance=gVariance, gByLocVar=gByLocVar, gByYearVar=gByYearVar, fracGxEAdd=fracGxEAdd, randLoc=randLoc, locCov=locCov)
    return(data)
  }
  with(sEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, loc.func)
      sfStop()
    }else{
      sims <- lapply(sims, loc.func)
    }
  })
}
