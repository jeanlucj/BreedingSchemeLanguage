#'Random mate
#'
#'@param simEnv an environment that BSL statements operate on
#'@param nProgeny the number of progenies
#'@param equalContribution if T all individuals used the same number of times as parents, if F individuals chosen at random to be parents
#'@param popID population ID to be crossed (default: the latest population)
#'@param popID2 population ID to be crossed with popID to make hybrids
#'
#'@return sequence information of progenies and the all information created before (list)
#'
#'@export
cross <- function(simEnv, nProgeny=100, equalContribution=F, popID=NULL, popID2=NULL){
  parent.env(simEnv) <- environment()
  cross.func <- function(data, nProgeny, equalContribution, popID, popID2){
    locPos <- data$mapData$map$Pos
    if(is.null(popID)){
      popID <- max(data$genoRec$popID)
    }
    tf <- data$genoRec$popID %in% popID
    GID.1 <- data$genoRec$GID[tf]
    nPar1 <- length(GID.1)
    geno <- data$geno[rep(GID.1*2, each=2) + rep(-1:0, nPar1), ]
    if (is.null(popID2)){
      if(equalContribution){
        geno <- randomMateAll(popSize=nProgeny, geno=geno, pos=locPos)
      }else{
        geno <- randomMate(popSize=nProgeny, geno=geno, pos=locPos)
      }
      pedigree <- matrix(GID.1[geno$pedigree], nrow=nProgeny)
      geno <- geno$progenies
    } else{ # Make pedigrees to mate two populations with each other
      tf <- data$genoRec$popID %in% popID2
      GID.2 <- data$genoRec$GID[tf]
      nPar2 <- length(GID.2)
      geno <- rbind(geno, data$geno[rep(GID.2*2, each=2) + rep(-1:0, nPar2), ])
      par1 <- c(rep(1:nPar1, nProgeny %/% nPar1), sample(nPar1, nProgeny %% nPar1))
      par2 <- nPar1 + c(rep(1:nPar2, nProgeny %/% nPar2), sample(nPar2, nProgeny %% nPar2))
      parents <- cbind(sample(par1), sample(par2))
      geno <- makeProgenies(parents, geno, locPos)
      pedigree <- cbind(GID.1[parents[,1]], GID.2[parents[,2]-nPar1])
    }
    GID <- max(data$genoRec$GID) + 1:nProgeny
    popID <- rep(max(data$genoRec$popID) + 1, nProgeny)
    hasGeno <- rep(FALSE, nProgeny)
    gValue <- calcGenotypicValue(geno=geno, mapData=data$mapData)
    addRec <- data.frame(GID=GID, pedigree=pedigree, popID=popID, basePopID=popID, hasGeno=hasGeno, gValue=gValue)
    data$genoRec <- rbind(data$genoRec, addRec)
    data$geno <- rbind(data$geno, geno)
    return(data)
  }
  with(simEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, cross.func, nProgeny=nProgeny, equalContribution=equalContribution, popID=popID, popID2=popID2)
      sfStop()
    } else{
      sims <- lapply(sims, cross.func, nProgeny=nProgeny, equalContribution=equalContribution, popID=popID, popID2=popID2)
    }
  })
}
