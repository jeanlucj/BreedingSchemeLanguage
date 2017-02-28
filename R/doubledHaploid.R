#'Doubled haploids
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param nProgeny the number of progeny
#'@param popID population ID to be devided by meiosis and doubled (default: the latest population)
#'
#'@return sequence information of progenies and the all information created before (list)
#'
#'@export
doubledHaploid <- function(sEnv=simEnv, nProgeny=100, popID=NULL){
  parent.env(sEnv) <- environment()
  doubledHaploid.func <- function(bsl, nProgeny, popID){
    locPos <- bsl$mapData$map$Pos
    if(is.null(popID)){
      popID <- max(bsl$genoRec$popID)
    }
    tf <- bsl$genoRec$popID %in% popID
    GIDpar <- bsl$genoRec$GID[tf]
    nPar <- length(GIDpar)
    geno <- bsl$geno[rep(GIDpar*2, each=2) + rep(-1:0, nPar),]
    geno <- makeDHs(popSize=nProgeny, geno=geno, pos=locPos)
    pedigree <- cbind(matrix(GIDpar[geno$pedigree], nProgeny), -1)
    geno <- geno$progenies
    bsl <- addProgenyData(bsl, geno, pedigree)
    if (exists("totalCost", bsl)) bsl$totalCost <- bsl$totalCost + nProgeny * bsl$costs$doubHapCost
    return(bsl)
  }
  with(sEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, doubledHaploid.func, nProgeny=nProgeny, popID=popID)
      sfStop()
    } else{
      sims <- lapply(sims, doubledHaploid.func, nProgeny=nProgeny, popID=popID)
    }
  })
}
