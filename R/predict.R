#'Genomic prediction
#'
#'@param simEnv an environment that BSL statements operate on
#'@param popID population ID to be predicted (default: the latest population)
#'@param trainingPopID population ID to be used for training a prediction model (default: all populations with phenotype data)
#'
#'@return predicted values and the all information created before (list)
#'
#'@export
predictBreedVal <- function(simEnv, popID = NULL, trainingPopID = NULL){
  parent.env(simEnv) <- environment()
  predict.func <- function(data, popID, trainingPopID){
    breedingData <- data$breedingData
    trainCandidates <- breedingData$hasGeno & breedingData$GID %in% breedingData$phenoGID
    if(is.null(trainingPopID)){
      GID.train <- sort(breedingData$GID[trainCandidates])
    }else{
      tf <- breedingData$popID %in% trainingPopID
      GID.train <- breedingData$GID[tf & trainCandidates]
    }
    mrkPos <- data$mapData$markerPos
    M <- (breedingData$geno[GID.train*2 - 1, mrkPos] + breedingData$geno[GID.train*2, mrkPos]) / 2
    R <- breedingData$error
    phenoGID <- breedingData$phenoGID
    y <- breedingData$pValue
    
    hetErr <- any(R != R[1])
    if (hetErr) sqrt.R <- sqrt(R)
    mt1ObsPerGID <- sum(phenoGID %in% GID.train) > length(GID.train)
    if (mt1ObsPerGID) factPhenGID <- as.factor(phenoGID)

    switch(
      hetErr + 2*mt1ObsPerGID + 1, 
      { # Homogeneous error, one obs per GID
        Z <- M
        X <- rep(1, length(y))
      },
      { # Heterogeneous error, one obs per GID
        y <- y / sqrt.R
        Z <- M / sqrt.R
        X <- 1 / sqrt.R # OK if X is a vector
      },
      { # Homogeneous error, more than one obs per GID
        w <- c(sqrt(table(phenoGID)))
        y <- c(tapply(y, factPhenGID, sum)) / w
        Z <- w * M
        X <- w
      },
      { # Heterogeneous error, more than one obs per GID
        w <- c(sqrt(tapply(1 / R, factPhenGID, sum)))
        y <- c(tapply(y / R, factPhenGID, sum)) / w
        Z <- w * M
        X <- c(tapply(1 / R, factPhenGID, sum)) / w
      }
    )
    
    mso <- tryCatch(mixed.solve(y, Z, X=X), error=function(e){
      print(e)
      print(Sys.time())
      save(data, y, Z, X, file="errorInMixedSolve.RData")
    })

    # Figure out who to predict
    if (is.null(popID)) popID <- max(breedingData$popID)
    tf <- breedingData$popID %in% popID
    GID.pred <- setdiff(breedingData$GID[tf], GID.train)
    predGID <- c(GID.train, GID.pred)
    M <- rbind(M, (breedingData$geno[GID.pred*2 - 1, mrkPos] + breedingData$geno[GID.pred*2, mrkPos]) / 2)
    predict <- M %*% mso$u

    if(is.null(breedingData$predNo)){
      predNo <- 1
    } else{
      predNo <- max(breedingData$predNo) + 1
    }
    breedingData$predict <- c(breedingData$predict, predict)
    breedingData$predGID <- c(breedingData$predGID, predGID)
    breedingData$predNo <- c(breedingData$predNo, rep(predNo, length(predGID)))
    
    selCriterion <- list(popID = popID, criterion = "pred")
    data$breedingData <- breedingData
    data$selCriterion <- selCriterion
    return(data)
  }#END predict.func
  
  with(simEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, predict.func, popID = popID, trainingPopID = trainingPopID)
      sfStop()
    }else{
      sims <- lapply(sims, predict.func, popID = popID, trainingPopID = trainingPopID)
    }
  })
}
