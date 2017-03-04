#' create map and QTL effects
#'
#' @param map map information (Chromosome and Position)
#' @param nLoci the number of markers and QTL
#' @param nMarkers the number of markers, which is used especially for genomic selection
#' @param nQTL the number of QTLs controlling the target trait
#' @param propDomi the probability of dominant QTL among the all QTL
#' @param interactionMean the expected number of epistatic loci for each effect
#' @param qtlInfo possible to give the function all of the qtl information
#'
#' @return map data including which loci are primary QTL and which are modifying loci, which have dominance effects, the effect sizes
#'
makeMap <- function(map, nLoci, nMarkers, nQTL, propDomi, interactionMean, qtlInfo=NULL){
  if (!is.null(qtlInfo)){
    actionType <- qtlInfo$actionType
    effectID <- qtlInfo$effectID
    effectivePos <- qtlInfo$effectivePos
    effects <- qtlInfo$effects
  } else {
    nEffectiveLoci <- 1 + rpois(n=nQTL, lambda=interactionMean)
    effectivePos <- sample(1:nLoci, sum(nEffectiveLoci))
    actionType <- rbinom(sum(nEffectiveLoci), 1, propDomi)
    effectID <- rep(1:nQTL, times=nEffectiveLoci)
    effects <- matrix(rnorm(nQTL), nQTL) # Will be modified by defineVariances if called
    if (nLoci - length(effectivePos) < nMarkers) stop("Number of markers is not enough!")
  }
  
  mrkPos <- sort(sample((1:nLoci)[-effectivePos], nMarkers, replace=F))
  return(list(map=map, markerPos=mrkPos, effectID=effectID, effectivePos=effectivePos, actionType=actionType, effects=effects))
}
