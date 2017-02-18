#'Doubled haploids
#'
#'@param simEnv an environment that BSL statements operate on
#'@param nProgeny the number of progeny
#'@param popID population ID to be devided by meiosis and doubled (default: the latest population)
#'
#'@return sequence information of progenies and the all information created before (list)
#'
#'@export
doubledHaploid <- function(simEnv, nProgeny=100, popID=NULL){
  parent.env(simEnv) <- environment()
  doubledHaploid.func <- function(data, nProgeny, popID){
    locPos <- data$mapData$map$Pos
    if(is.null(popID)){
      popID <- max(data$genoRec$popID)
    }
    tf <- data$genoRec$popID %in% popID
    GIDpar <- data$genoRec$GID[tf]
    nPar <- length(GIDpar)
    geno <- data$geno[rep(GIDpar*2, each=2) + rep(-1:0, nPar),]
    geno <- makeDHs(popSize=nProgeny, geno=geno, pos=locPos)
    pedigree <- matrix(GIDpar[geno$pedigree], nProgeny)
    geno <- geno$progenies
    gValue <- calcGenotypicValue(geno=geno, mapData=data$mapData)
    GID <- max(data$genoRec$GID) + 1:nProgeny
    popID <- max(data$genoRec$popID) + 1
    gValue <- calcGenotypicValue(geno=geno, mapData=data$mapData)
    addRec <- data.frame(GID, pedigree, popID, basePopID=popID, hasGeno=FALSE, gValue)
    colnames(addRec) <- colnames(data$genoRec)
    data$genoRec <- rbind(data$genoRec, addRec)
    data$geno <- rbind(data$geno, geno)
    return(data)
  }
  with(simEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, doubledHaploid.func, nProgeny=nProgeny, popID=popID)
      sfStop()
    } else{
      sims <- lapply(sims, doubledHaploid.func, nProgeny=nProgeny, popID=popID)
    }
  })
}
