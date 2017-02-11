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
    makeGamete(geno[par * 2 + 0:1, ], pos)
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
  parent <- rep(1:nPar, nRep)
  parent <- c(parent, sample(1:nPar, rem))
  progenies <- t(sapply(parent, function(par) DH(geno[c(parent[par] * 2 - 1, parent[par] * 2), ], pos)))
  return(list(progenies = progenies, pedigree = cbind(parent, parent)))
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
  parent <- rep(1:nPar, nRep)
  parent <- c(parent, sample(1:nPar, rem))
  progenies <- makeProgenies(cbind(parent, parent), geno, pos)
  return(list(progenies = progenies, pedigree = cbind(parent, parent)))
}

#'randomMate
#'
#'@param popSize population size
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
randomMate <- function(popSize, geno, pos){
  parents <- t(sapply(rep(nrow(geno) / 2, popSize), sample, size=2))
  progenies <- makeProgenies(parents, geno, pos)
  return(list(progenies = progenies, pedigree = parents))
}

#'randomMateAll
#'
#'@param popSize population size
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
randomMateAll <- function(popSize, geno, pos){
  nInd <- nrow(geno) / 2
  parent1 <- rep(1:nInd, length.out=popSize)
  parent2 <- sapply(parent1, function(par) sample((1:nInd)[-parent1[par]], size = 1))
  parents <- cbind(parent1, parent2)
  progenies <- makeProgenies(parents, geno, pos)
  return(list(progenies = progenies, pedigree = parents))
}
