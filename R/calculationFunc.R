#'calcGenotypicValue
#'
#'@param geno matrix of haplotypes
#'@param mapData map data
#'
# Calculate genotypic value one QTL at a time for all individuals
calcGenotypicValue <- function(geno, mapData){
  # Calculate the contribution to genotypic value of one indiviudal of one QTL
  gv1indThisQ <- function(ind){ # ind is zero base
    geno1pos <- geno[ind*2 + 1:2, posThisQ, drop=F]
    if (!is.null(mapData$domModel)){ # 1 dominant over -1; degree in actType
      minMaxGeno <- apply(geno1pos, 2, range)
      coef <- sapply(1:length(actType), function(i) c(1-actType[i], actType[i]) %*% minMaxGeno[,i])
    } else{ # Standard model
      coef <- ifelse(actType == 0, (geno1pos[1, ] + geno1pos[2, ])/2, -(geno1pos[1, ] * geno1pos[2, ]))
    }
    return(effect * prod(coef))
  }

  nInd <- nrow(geno) / 2
  genoVal <- numeric(nInd)
  for(i in 1:max(mapData$effectID)){
    posThisQ <- mapData$effectivePos[mapData$effectID == i]
    actType <- mapData$actionType[mapData$effectID == i]
    effect <- mapData$effects[i, 1]
    genoVal <- genoVal + sapply(0:(nInd - 1), gv1indThisQ)
  }
  return(genoVal)
}

#'calcPhenotypicValue
#'
#'@param gv genotypic values
#'@param errorVar error variance
#'@param H2 a broad heritability
#'
calcPhenotypicValue <- function(gv, errorVar = NULL, H2 = NULL){
  if((!is.null(errorVar) & !is.null(H2)) | (is.null(errorVar) & is.null(H2))){
    stop("I cannot make phenotypic value!")
  }else{
    if(!is.null(errorVar)){
      pv <- gv + rnorm(length(gv), 0, sqrt(errorVar))
    }else{
      varG <- var(gv)
      errorVar <- varG * (1 - H2) / H2
      pv <- gv + rnorm(length(gv), 0, sqrt(errorVar))
    }
  }
  return(pv)
}
