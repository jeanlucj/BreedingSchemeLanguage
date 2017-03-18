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
  crossOver <- rec >= runif(length(rec))
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
#'@param popSize population size
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
#'@param popSize population size
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
# Randomly mate with no selfing
randomMate <- function(popSize, geno, pos){
  parents <- t(sapply(rep(nrow(geno) / 2, popSize), sample, size=2))
  progenies <- makeProgenies(parents, geno, pos)
  return(list(progenies = progenies, pedigree = parents))
}

# Randomly mate but all parents have to be used equally.
# It's trickier than it seems
#'randomMateAll
#'
#'@param popSize population size
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
randomMateAll <- function(popSize, geno, pos){
  equalAndRand <- function(popSize, nPar){
    parents <- matrix(sample(c(rep(1:nPar, 2*popSize %/% nPar), sample(nPar, 2*popSize %% nPar))), popSize)
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
