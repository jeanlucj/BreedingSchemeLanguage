#'makeGamete
#'
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
makeGamete <- function(geno, pos){
  btwLocDist <- diff(pos)
  rec <- (1 - exp(-2 * btwLocDist / 100)) / 2
  rec[rec < 0] <- 0.5
  rec <- c(0.5, rec)
  crossOver <- rec >= stats::runif(length(rec))
  selectHaplo <- cumsum(crossOver) %% 2
  return(ifelse(selectHaplo, geno[1,], geno[2,]))
}

#'makeProgeny
#'
#'@param genoPat matrix of paternal haplotype
#'@param genoMat matrix of maternal haplotype
#'@param pos position of markers/QTLs
#'
makeProgeny <- function(genoMat, genoPat, pos){
  return(rbind(makeGamete(genoMat, pos), makeGamete(genoPat, pos)))
}

#'makeProgenies
#'
#'@param parents ID of haplotypes that the parents harbor
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
makeProgenies <- function(parents, geno, pos){
  gameteOnePar <- function(par){
    makeGamete(geno[par * 2 + -1:0, ], pos)
  }
  return(t(sapply(c(t(parents)), gameteOnePar)))
}

#'DH
#'
#'@param genoParent matrix of haplotypes
#'@param pos position of markers/QTLs
#'
DH <- function(genoParent, pos){
  gamete <- makeGamete(genoParent, pos)
  progeny <- rbind(gamete, gamete)
  return(progeny)
}

#'makeDHs
#'
#'@param popSize the number of DH individuals to return
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
makeDHs <- function(popSize, geno, pos){
  nPar <- nrow(geno) / 2
  nRep <- popSize %/% nPar
  rem <- popSize %% nPar
  parents <- c(rep(1:nPar, nRep), sample(1:nPar, rem))
  progenies <- t(sapply(parents, function(par) makeGamete(geno[par*2 + -1:0, ], pos)))
  progenies <- rbind(progenies, progenies)[rep(c(0, popSize), popSize) + rep(1:popSize, each=2), ]
  return(list(progenies = progenies, pedigree = cbind(parents, parents)))
}

#'makeSelfs
#'
#'@param popSize the number of selfed individuals to return
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
makeSelfs <- function(popSize, geno, pos){
  nPar <- nrow(geno) / 2
  nRep <- popSize %/% nPar
  rem <- popSize %% nPar
  parents <- c(rep(1:nPar, nRep), sample(1:nPar, rem))
  parents <- cbind(parents, parents)
  progenies <- makeProgenies(parents, geno, pos)
  return(list(progenies=progenies, pedigree=parents))
}

#'randomMate
#'
#'@param popSize the number of progeny to return
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
# Randomly mate. Selfing is possible.
randomMate <- function(popSize, geno, pos){
  parents <- t(sapply(rep(nrow(geno) / 2, popSize), sample, size=2, replace=T))
  progenies <- makeProgenies(parents, geno, pos)
  return(list(progenies = progenies, pedigree = parents))
}

#'pedigreeMate
#'
#'@param parents two- or three-column matrix: first two columns index the parents you want to cross, if there is a third column, it is the number of progeny from that pair of parents
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
# Make progeny based on a pedigree.
pedigreeMate <- function(parents, geno, pos){
  if (ncol(parents) > 2){
    test <- apply(parents, 1, function(vec) rep.int(vec[1:2], vec[3]))
    test <- matrix(unlist(test), ncol=2, byrow=T)
  }
  progenies <- makeProgenies(parents, geno, pos)
  return(list(progenies = progenies, pedigree = parents))
}

# Randomly mate but all parents have to be used equally and no selfing is allowed.
# It's trickier than it seems
#'randomMateAll
#'
#'@param popSize the number of progeny to return
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
randomMateAll <- function(popSize, geno, pos){
  equalAndRand <- function(popSize, nPar){
    parents <- matrix(sample(c(rep(1:nPar, (2*popSize) %/% nPar), sample(nPar, (2*popSize) %% nPar))), popSize)
    noSelfs <- function(parRow){
      if (parents[parRow, 1] == parents[parRow, 2]){
        par <- parents[parRow, 1]
        swapCan <- apply(parents, 1, function(can) sum(can == par))
        swapRow <- sample(which(swapCan == 0), 1)
        parents[parRow, 1] <<- parents[swapRow, 1]
        parents[swapRow, 1] <<- par 
      }
    }
    dummy <- sapply(1:popSize, noSelfs)
    return(parents)
  }
  
  nPar <- nrow(geno) / 2
  parents <- equalAndRand(popSize, nPar)
  progenies <- makeProgenies(parents, geno, pos)
  return(list(progenies = progenies, pedigree = parents))
}

# Randomly mate but use parents equally and do not allow matings between full- or half-sibs
#'randomMateNoFam
#'
#'@param popSize the number of progeny to return
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'@param genoRec pedigree records of the individuals represented by geno
#'
randomMateNoFam <- function(popSize, geno, pos, genoRec){
  equalAndNoFam <- function(popSize, nPar, genoRec){
    parents <- matrix(sample(c(rep(1:nPar, (2*popSize) %/% nPar), sample(nPar, (2*popSize) %% nPar))), popSize)
    noFamPairs <- function(parRow){
      related <- function(par1, par2){
        any(genoRec[par1, 2:3] %in% genoRec[par2, 2:3]) # Columns 2 and 3 have the parent IDs
      }
      okSwap <- function(pair){
        cond1 <- sapply(parents[,1], related, par2=pair[2])
        cond2 <- sapply(parents[,2], related, par2=pair[1])
        return(which(!cond1 & !cond2))
      }
      if (related(parents[parRow, 1], parents[parRow, 2])){
        par <- parents[parRow, 1]
        swapCan <- okSwap(parents[parRow,])
        swapRow <- sample(swapCan, 1)
        parents[parRow, 1] <<- parents[swapRow, 1]
        parents[swapRow, 1] <<- par 
      }
    }
    dummy <- sapply(1:popSize, noFamPairs)
    return(parents)
  }
  
  nPar <- nrow(geno) / 2
  parents <- equalAndNoFam(popSize, nPar, genoRec)
  progenies <- makeProgenies(parents, geno, pos)
  return(list(progenies = progenies, pedigree = parents))
}

#'predGameteMeanVar
#'
#'@param geno 2 x nLoc matrix of haplotypes, coded -1 and 1
#'@param pos position of markers/QTLs
#'@param locEff effects of each locus
#'
predGameteMeanVar <- function(geno, pos, locEff){
  expMean <- colMeans(geno) %*% locEff
  polymorphic <- apply(geno, 2, stats::sd) > 0
  nPoly <- sum(polymorphic)
  pos <- pos[polymorphic]; locEff <- locEff[polymorphic]
  btwLocDist <- diff(pos)
  chrBrk <- which(btwLocDist < 0) + 1
  sameChr <- matrix(1, nPoly, nPoly)
  for (b in chrBrk){
    sameChr[1:(b-1), b:nPoly] <- sameChr[b:nPoly, 1:(b-1)] <- 0
  }
  gamDiseq <- exp(-2 * as.matrix(stats::dist(btwLocDist)) / 100)
  gamDiseq <- gamDiseq * sameChr
  signSqrt <- sqrt(abs(locEff))*sign(locEff)
  expVar <- t(signSqrt) %*% gamDiseq %*% signSqrt
  return(c(expMean=expMean, expVar=expVar))
}
