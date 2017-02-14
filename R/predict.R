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
    mt1ObsPerGID <- sum(phenoGID %in% GID.train) > length(GID.train)

    # Figure out who to predict
    if (is.null(popID)) popID <- max(breedingData$popID)
    tf <- breedingData$popID %in% popID
    GID.pred <- setdiff(breedingData$GID[tf], GID.train)
    M <- rbind(M, (breedingData$geno[GID.pred*2 - 1, mrkPos] + breedingData$geno[GID.pred*2, mrkPos]) / 2)
    predGID <- c(GID.train, GID.pred)
    
    K <- A.mat(M)
    rownames(K) <- colnames(K) <- predGID
    keep <- breedingData$phenoGID %in% predGID
    kbDat <- data.frame(breedingData$pValue[keep], breedingData$phenoGID[keep])
    colnames(kbDat) <- c("pheno", "GID")
    kbo <- kin.blup(kbDat, geno="GID", pheno="pheno", K=K, reduce=mt1ObsPerGID, R=R)
    predict <- kbo$g[as.character(predGID)]
    
    if(is.null(breedingData$predict)){
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
