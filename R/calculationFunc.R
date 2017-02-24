#'calcGenotypicValue
#'
#'@param geno matrix of haplotypes
#'@param mapData map data
#'
# Calculate the genotypic values for individuals with geno at all locations
calcGenotypicValue <- function(geno, mapData){
  # Calculate genotypic value one QTL at a time for all individuals
  gv1indThisQ <- function(ind){ # ind is zero base
    geno1pos <- geno[ind*2 + 1:2, posThisQ, drop=F]
    if (!is.null(mapData$domModel)){ # 1 dominant over -1; degree in actType
      minMaxGeno <- apply(geno1pos, 2, range)
      coef <- sapply(1:length(actType), function(i) c(1-actType[i], actType[i]) %*% minMaxGeno[,i])
    } else{ # Standard model (-1,-1 same as 1,1 and opposite to -1,1 or 1,-1)
      coef <- ifelse(actType == 0, (geno1pos[1, ] + geno1pos[2, ])/2, -(geno1pos[1, ] * geno1pos[2, ]))
    }
    return(effect * prod(coef))
  }#END gv1indThisQ

  nInd <- nrow(geno) / 2
  genoVal <- matrix(0, ncol(mapData$effects), nInd)
  for(i in 1:max(mapData$effectID)){
    posThisQ <- mapData$effectivePos[mapData$effectID == i]
    actType <- mapData$actionType[mapData$effectID == i]
    effect <- mapData$effects[i, ]
    genoVal <- genoVal + sapply(0:(nInd - 1), gv1indThisQ)
  }
  return(t(genoVal))
}

#'calcPhenotypicValue
#'
#'@param gv genotypic values
#'@param errorVar error variance
#'
# If gv is a multitrait matrix, errorVar is equal across all traits
calcPhenotypicValue <- function(gv, errorVar){
  pv <- gv + rnorm(length(gv), 0, sqrt(errorVar))
  return(pv)
}
