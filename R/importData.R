#' Transform a data.frame with a hapmap data in it into a marker dosage and map list
#' 
#' @param hm The data.frame that has HapMap data in it
#' 
#' @return list with markers and map objects
#' 
# Assume marker names in column 1 and genotypes (e.g. "A/G") in column 2
# Assume chromosome and position are in columns 3 and 4 of the hapmap
# and that position is in cM: WARNING cM is probably not standard hapmap
phasedHapMap2mat <- function(hm){
  map <- data.frame(Chr=hm[[3]], Pos=hm[[4]])
  # See if user specified QTL information, and if so pull it in
  # That information would be coded in columns 7, 8, and 9
  # col7 = effectID, col8 = actionType, col9 = effect
  if (any(!is.na(c(hm[[7]], hm[[8]], hm[[9]])))){
    qtl <- which(!is.na(hm[[8]]))
    actionType <- hm[[8]][qtl]
    effects <- hm[[9]][qtl]
    effects <- matrix(effects[!is.na(effects)], ncol=1)
    if (all(!is.na(hm[[7]][qtl]))){
      effectID <- as.integer(hm[[7]][qtl])
    } else {
      effectID <- 1:length(qtl)
    }
    qtlInfo <- list(effectID=effectID, effectivePos=qtl, actionType=actionType, effects=effects)
  } else qtlInfo <- NULL
  hm <- apply(hm[c(2, 12:ncol(hm))], 2, as.character)
  # Different functions if hapmap with 1 or 2 character codes
  if (nchar(hm[1, 2]) == 1){ 
    hapMap2num <- function(vec){ # Convert to numeric if one character codes
      hets <- c("W","R","Y","K","M","S") # NOTE: hets become missing here
      missings <- c("-", "0","N")
      alleles <- unlist(strsplit(vec[1], "/"))
      vec <- vec[-1]
      vec[vec %in% c(hets, missings)] <- NA
      vec[alleles[1]==vec] <- 1
      vec[alleles[2]==vec] <- 0
      return(suppressWarnings(as.numeric(vec)))
    }
  } else{ # Convert to numeric if two character codes
    hapMap2num <- function(vec){
      missings <- c("-", "0","N")
      alleles <- unlist(strsplit(vec[1], "/"))
      vec <- vec[-1]
      codes <- c(alleles, missings)
      numCodes <- c(1, 0, rep(NA, length(missings)))
      gam1 <- sapply(substr(vec, 1, 1), function(code) numCodes[code == codes])
      gam2 <- sapply(substr(vec, 2, 2), function(code) numCodes[code == codes])
      return(c(gam1, gam2))
    }
  }
  res <- apply(hm, 1, hapMap2num) 
  return(list(markers=res, map=map, qtlInfo=qtlInfo))
}

#' Generate a data.frame with a hapmap data in it to test phasedHapMap2mat
#' 
#' @param nInd The number of individuals with marker data
#' @param nMrk The number of markers
#' @param nChr The number of chromosomes the species has
#' @param lenChr The length of the chromosomes (all equal) in cM
#' @param maf The desired minor allele frequency of each marker
#' @param nCharCode Whether the genotype codes should be one or two characters
#' @param nQTL If you want to put QTL information in the hapmap
#' @param propDomi proportion of QTL that act dominantly
#' @param interactionMean average number of loci interacting with QTL: Poisson
#' @param varEffects variance of the QTL effects
#' 
#' @return A data.frame with hapmap data in it
#' 
# Quick function to create a HapMap data.frame for testing
#
simHapMap <- function(nInd=200, nMrk=1050, nChr=7, lenChr=150, maf=runif(nMrk), nCharCode=2, nQTL=NULL, propDomi=NULL, interactionMean=NULL, varEffects=NULL){
  nucl <- c("A", "C", "G", "T")
  gt <- replicate(nMrk, sample(4, 2))
  alleles <- paste(nucl[gt[1,]], nucl[gt[2,]], sep="/")
  chr <- sort(sample(nChr, nMrk, replace=T))
  pos <- round(lenChr * runif(nMrk), 2)
  pos <- pos[order(chr, pos)]
  hapmap <- data.frame(locus=paste("locus", 1:nMrk, sep=""), alleles=alleles, chrom=chr, pos=pos)
  if (!(nCharCode %in% 1:2)) nCharCode <- 2
  nInd <- nInd * 2 / nCharCode
  gam <- matrix(rbinom(nInd*nMrk, 1, 1 - maf) + 1, nMrk, nInd)
  gam1 <- apply(gam, 2, function(gamVec) nucl[gt[cbind(gamVec, 1:nMrk)]])
  if (nCharCode == 2){
    gam <- matrix(rbinom(nInd*nMrk, 1, 1 - maf) + 1, nMrk, nInd)
    gam2 <- apply(gam, 2, function(gamVec) nucl[gt[cbind(gamVec, 1:nMrk)]])
  } else gam2 <- NULL
  hm <- matrix(paste(gam1, gam2, sep=""), nMrk, nInd)
  if (!is.null(nQTL)){
    nEffectiveLoci <- 1 + rpois(n = nQTL, lambda = interactionMean)
    effectivePos <- sample(1:nMrk, sum(nEffectiveLoci))
    actionType <- rbinom(sum(nEffectiveLoci), 1, propDomi)
    effectID <- rep(1:nQTL, times=nEffectiveLoci)
    effects <- as.matrix(rnorm(nQTL, 0, sqrt(varEffects)), nQTL)
    col7 <- col8 <- col9 <- rep(NA, nMrk)
    col7[effectivePos] <- effectID
    col8[effectivePos] <- actionType
    idx <- 1
    for (q in 1:nQTL){
      col9[effectivePos[idx]] <- effects[q]
      idx <- idx + nEffectiveLoci[q]
    }
    hapmap <- cbind(hapmap, matrix(NA, nMrk, 2), effectID=col7, actionType=col8, effects=col9, matrix(NA, nMrk, 2), hm)
  } else hapmap <- cbind(hapmap, matrix(NA, nMrk, 7), hm)
  colnames(hapmap)[12:(nInd+11)] <- paste("ind", 1:nInd, sep="")
  return(hapmap)
}
