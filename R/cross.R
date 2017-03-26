#' Cross with random mating, or equal contributions,
#'or randomly between two populations
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param nProgeny the number of progenies
#'@param equalContribution if T all individuals used the same number of times as parents, if F individuals chosen at random to be parents
#'@param popID population ID to be crossed. Default: the last population
#'@param popID2 population ID to be crossed with popID to make hybrids. Default: second population not used.
#'
#'@return sequence information of progenies and the all information created before (list)
#'
#'@export
cross <- function(sEnv=simEnv, nProgeny=100, equalContribution=F, popID=NULL, popID2=NULL){
  parent.env(sEnv) <- environment()
  cross.func <- function(bsl, nProgeny, equalContribution, popID, popID2){
    locPos <- bsl$mapData$map$Pos
    if(is.null(popID)){
      popID <- max(bsl$genoRec$popID)
    }
    tf <- bsl$genoRec$popID %in% popID
    GID.1 <- bsl$genoRec$GID[tf]
    nPar1 <- length(GID.1)
    geno <- bsl$geno[rep(GID.1*2, each=2) + rep(-1:0, nPar1), ]
    if (is.null(popID2)){
      if(equalContribution){
        geno <- randomMateAll(popSize=nProgeny, geno=geno, pos=locPos)
      }else{
        geno <- randomMate(popSize=nProgeny, geno=geno, pos=locPos)
      }
      pedigree <- cbind(matrix(GID.1[geno$pedigree], nrow=nProgeny), 0)
      geno <- geno$progenies
    } else{ # Make pedigrees to mate two populations with each other
      tf <- bsl$genoRec$popID %in% popID2
      GID.2 <- bsl$genoRec$GID[tf]
      nPar2 <- length(GID.2)
      geno <- rbind(geno, bsl$geno[rep(GID.2*2, each=2) + rep(-1:0, nPar2), ])
      par1 <- sample(c(rep(1:nPar1, nProgeny %/% nPar1), sample(nPar1, nProgeny %% nPar1)))
      par2 <- nPar1 + sample(c(rep(1:nPar2, nProgeny %/% nPar2), sample(nPar2, nProgeny %% nPar2)))
      parents <- cbind(sample(par1), sample(par2))
      geno <- makeProgenies(parents, geno, locPos)
      pedigree <- cbind(GID.1[parents[,1]], GID.2[parents[,2]-nPar1], 0)
    }
    bsl <- addProgenyData(bsl, geno, pedigree)
    if (exists("totalCost", bsl)) bsl$totalCost <- bsl$totalCost + nProgeny * bsl$costs$crossCost
    return(bsl)
  }#END cross.func

  with(sEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, cross.func, nProgeny=nProgeny, equalContribution=equalContribution, popID=popID, popID2=popID2)
      sfStop()
    } else{
      sims <- lapply(sims, cross.func, nProgeny=nProgeny, equalContribution=equalContribution, popID=popID, popID2=popID2)
    }
  })
}

#' Add progeny information to data after cross, doubledHaploid, or selfFertilize
#'
#'@param bsl the list that has all the objects for one simulation
#'@param geno the genotypes of the progeny
#'@param pedigree the three-column pedigree of the progeny (last col: DH, outbred, self)
#'
#'@return data with progeny information added
#'
addProgenyData <- function(bsl, geno, pedigree){
  # Add on to genetic values
  gValue <- calcGenotypicValue(geno=geno, mapData=bsl$mapData)
  bsl$gValue <- rbind(bsl$gValue, gValue)
  # Add on to the QTL relationship matrix
  nProgeny <- nrow(geno) / 2
  M <- geno[1:nProgeny*2 - 1, bsl$mapData$effectivePos] + geno[1:nProgeny*2, bsl$mapData$effectivePos]
  nAdd <- ncol(bsl$yearEffects)
  if (nAdd == 0){ # The user is making croses without ever having phenotyped
    bsl$locEffects <- matrix(0, nrow=nrow(bsl$locEffects) + nProgeny, ncol=0)
    bsl$locEffectsI <- matrix(0, nrow=nrow(bsl$locEffectsI) + nProgeny, ncol=0)
    bsl$yearEffects <- matrix(0, nrow=nrow(bsl$yearEffects) + nProgeny, ncol=0)
    bsl$yearEffectsI <- matrix(0, nrow=nrow(bsl$yearEffectsI) + nProgeny, ncol=0)
  } else{
    vp <- bsl$varParms$gByYearVar * bsl$varParms$fracGxEAdd
    toAdd <- M %*% bsl$gByYqtl
    toAdd <- sapply(1:nAdd, function(i) toAdd[,i] * bsl$yearScale[i])
    bsl$yearEffects <- rbind(bsl$yearEffects, toAdd)
    vp <- bsl$varParms$gByYearVar * (1 - bsl$varParms$fracGxEAdd)
    toAdd <- matrix(rnorm(nProgeny * nAdd, sd=sqrt(vp)), nProgeny)
    bsl$yearEffectsI <- rbind(bsl$yearEffectsI, toAdd)
  }
  nAdd <- ncol(bsl$locEffects)
  if (bsl$varParms$randLoc & nAdd > 0){
    vp <- bsl$varParms$gByLocVar * bsl$varParms$fracGxEAdd
    toAdd <- M %*% bsl$gByLqtl
    toAdd <- sapply(1:nAdd, function(i) toAdd[,i] * bsl$locScale[i])
    bsl$locEffects <- rbind(bsl$locEffects, toAdd)
    vp <- bsl$varParms$gByLocVar * (1 - bsl$varParms$fracGxEAdd)
    toAdd <- matrix(rnorm(nProgeny * nAdd, sd=sqrt(vp)), nProgeny)
    bsl$locEffectsI <- rbind(bsl$locEffectsI, toAdd)
  }
  # Add on to the genotypic records
  GID <- max(bsl$genoRec$GID) + 1:nProgeny
  popID <- rep(max(bsl$genoRec$popID) + 1, nProgeny)
  hasGeno <- rep(FALSE, nProgeny)
  addRec <- data.frame(GID=GID, pedigree=pedigree, popID=popID, basePopID=popID, hasGeno=hasGeno)
  colnames(addRec) <- colnames(bsl$genoRec)
  bsl$genoRec <- rbind(bsl$genoRec, addRec)
  # Add on to the genotypes
  bsl$geno <- rbind(bsl$geno, geno)
  return(bsl)
}