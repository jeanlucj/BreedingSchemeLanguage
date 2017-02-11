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
    breedingData <- data$breedingData
    selCriterion <- data$selCriterion
    if(is.null(popID)){
      popID <- max(breedingData$popID)
    }
    tf <- breedingData$popID %in% popID
    GID.now <- breedingData$GID[tf]
    geno.now <- breedingData$geno[sort(c(GID.now * 2 - 1, GID.now * 2)), ]
    if (is.null(popID2)){
      if(equalContribution){
        geno.progeny <- randomMateAll(popSize=nProgeny, geno=geno.now, pos=locPos)$progenies
      }else{
        geno.progeny <- randomMate(popSize=nProgeny, geno=geno.now, pos=locPos)$progenies
      }
    } else{ # Make pedigrees to mate two populations with each other
      nPop1 <- sum(tf)
      tf <- breedingData$popID %in% popID2
      GID.2 <- breedingData$GID[tf]
      nPop2 <- sum(tf)
      geno.now <- rbind(geno.now, breedingData$geno[sort(c(GID.2 * 2 - 1, GID.2 * 2)), ])
      parents <- cbind(sample(rep(1:nPop1, length.out=nProgeny)), sample(rep(nPop1 + 1:nPop2, length.out=nProgeny)))
      geno.progeny <- makeProgenies(parents, geno.now, locPos)
    }
    gValue <- calcGenotypicValue(geno=geno.progeny, mapData=data$mapData)
    GID.progeny <- max(breedingData$GID) + 1:nProgeny
    breedingData$GID <- c(breedingData$GID, GID.progeny)
    popID.progeny <- rep(max(breedingData$popID) + 1, nProgeny)
    breedingData$popID <- c(breedingData$popID, popID.progeny)
    breedingData$popIDsel <- c(breedingData$popIDsel, popID.progeny)
    breedingData$geno <- rbind(breedingData$geno, geno.progeny)
    breedingData$hasGeno <- c(breedingData$hasGeno, rep(FALSE, nProgeny))
    breedingData$gValue <- c(breedingData$gValue, gValue)
    data$breedingData <- breedingData
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
