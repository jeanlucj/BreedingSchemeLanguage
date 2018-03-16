#' Cross with random mating, or equal contributions, or randomly between two populations
#'
#'@param sEnv the environment that BSL functions operate in. Default: "simEnv" so use that to avoid specifying when calling functions
#'@param nProgeny the number of progenies. Default: 100
#'@param equalContribution if T all individuals used the same number of times as parents, if F individuals chosen at random to be parents. Default: F
#'@param popID population ID to be crossed. Default: the last population
#'@param popID2 population ID to be crossed with popID to make hybrids. Default: second population not used
#'@param notWithinFam if TRUE, like equalContribution, all individuals are used the same number of times as parents and self-fertilization is not allowed. In addition, half- and full-sibs are not allowed to mate. Default: F
#'@param pedigree optional two- or three-column matrix: the first two columns are the GIDs that you want to cross, the third column is the number of progeny from that cross. NOTE: pedigree supercedes the nProgeny, equalContribution, popID, popID2, and notWithinFam parameters. You have to know what you are doing to use this parameter. Default: NULL
#'@param parms optional named list. Objects with those names will be created with the corresponding values. A way to pass values that are not predetermined by the script. Default: NULL
#'
#'@return modifies the list sims in environment sEnv by creating a progeny population as specified, with an incremented population number
#'
#'@export
cross <- function(sEnv=NULL, nProgeny=100, equalContribution=F, popID=NULL, popID2=NULL, notWithinFam=F, pedigree=NULL, parms=NULL){
  if(!is.null(parms)){
    for (n in 1:length(parms)){
      assign(names(parms)[n], parms[[n]])
    }
  }
  cross.func <- function(bsl, nProgeny, equalContribution, popID, popID2){
    locPos <- bsl$mapData$map$Pos
    if (!is.null(pedigree)){
      GID.1 <- sort(unique(c(pedigree[,1:2])))
      geno <- bsl$geno[sort(c(GID.1 * 2 - 1, GID.1 * 2)),]
      idx <- 1:length(GID.1)
      names(idx) <- GID.1
      parents <- matrix(idx[as.character(pedigree)], nrow(pedigree))
      if (ncol(pedigree) > 2) parents <- cbind(parents, pedigree[,3])
      geno <- pedigreeMate(parents=parents, geno=geno, pos=locPos)
      pedigree <- cbind(matrix(GID.1[geno$pedigree], ncol=2), 0)
      geno <- geno$progenies
    } else{
      if(is.null(popID)){
        popID <- max(bsl$genoRec$popID)
      }
      tf <- bsl$genoRec$popID %in% popID
      GID.1 <- bsl$genoRec$GID[tf]
      nPar1 <- length(GID.1)
      if (nPar1 == 0) stop(paste("There are no parents in population", popID))
      geno <- bsl$geno[rep(tf, each=2),]
      if (is.null(popID2)){
        if (equalContribution | notWithinFam){
          if (notWithinFam){ # notWithinFam fulfills the conditions of equalContribution
            geno <- randomMateNoFam(popSize=nProgeny, geno=geno, pos=locPos, bsl$genoRec[tf,])
          } else{
            geno <- randomMateAll(popSize=nProgeny, geno=geno, pos=locPos)
          }
        } else{
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
    }#END no pedigree was given
    bsl <- addProgenyData(bsl, geno, pedigree)
    if (exists("totalCost", bsl)) bsl$totalCost <- bsl$totalCost + nProgeny * bsl$costs$crossCost
    return(bsl)
  }#END cross.func
  
  if(is.null(sEnv)){
    if(exists("simEnv", .GlobalEnv)){
      sEnv <- get("simEnv", .GlobalEnv)
    } else{
      stop("No simulation environment was passed")
    }
  } 
  parent.env(sEnv) <- environment()
  with(sEnv, {
    if(nCore > 1){
      snowfall::sfInit(parallel=T, cpus=nCore)
      sims <- snowfall::sfLapply(sims, cross.func, nProgeny=nProgeny, equalContribution=equalContribution, popID=popID, popID2=popID2)
      snowfall::sfStop()
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
#' Normally DH=-1, outbred=0, self=nGen of selfing but that is not currenly properly updated
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
    toAdd <- matrix(stats::rnorm(nProgeny * nAdd, sd=sqrt(vp)), nProgeny)
    bsl$yearEffectsI <- rbind(bsl$yearEffectsI, toAdd)
  }
  nAdd <- ncol(bsl$locEffects)
  if (bsl$varParms$randLoc & nAdd > 0){
    vp <- bsl$varParms$gByLocVar * bsl$varParms$fracGxEAdd
    toAdd <- M %*% bsl$gByLqtl
    toAdd <- sapply(1:nAdd, function(i) toAdd[,i] * bsl$locScale[i])
    bsl$locEffects <- rbind(bsl$locEffects, toAdd)
    vp <- bsl$varParms$gByLocVar * (1 - bsl$varParms$fracGxEAdd)
    toAdd <- matrix(stats::rnorm(nProgeny * nAdd, sd=sqrt(vp)), nProgeny)
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