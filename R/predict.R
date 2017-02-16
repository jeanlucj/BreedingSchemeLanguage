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
    trainCandidates <- data$genoRec$hasGeno & data$genoRec$GID %in% data$phenoRec$phenoGID
    if(is.null(trainingPopID)){
      GID.train <- data$genoRec$GID[trainCandidates]
    }else{
      tf <- data$genoRec$popID %in% trainingPopID
      GID.train <- data$genoRec$GID[tf & trainCandidates]
    }
    mrkPos <- data$mapData$markerPos

    M <- (data$geno[GID.train*2 - 1, mrkPos] + data$geno[GID.train*2, mrkPos]) / 2
    R <- data$phenoRec$error
    phenoGID <- data$phenoRec$phenoGID
    y <- data$phenoRec$pValue
    mt1ObsPerGID <- sum(phenoGID %in% GID.train) > length(GID.train)

    # Figure out who to predict
    if (is.null(popID)) popID <- max(data$genoRec$popID)
    tf <- data$genoRec$popID %in% popID
    GID.pred <- setdiff(data$genoRec$GID[tf], GID.train)
    M <- rbind(M, (data$geno[GID.pred*2 - 1, mrkPos] + data$geno[GID.pred*2, mrkPos]) / 2)
    predGID <- c(GID.train, GID.pred)
    
    K <- A.mat(M)
    rownames(K) <- colnames(K) <- predGID
    keep <- data$phenoRec$phenoGID %in% predGID
    kbDat <- data.frame(pheno=data$phenoRec$pValue[keep], GID=data$phenoRec$phenoGID[keep])
    kbo <- kin.blup(kbDat, geno="GID", pheno="pheno", K=K, reduce=mt1ObsPerGID, R=R)
    predict <- kbo$g[as.character(predGID)]
    
    if(is.null(data$predRec)){
      predNo <- 1
    } else{
      predNo <- max(data$predRec$predNo) + 1
    }
    toAdd <- data.frame(predGID=predGID, predNo=predNo, predict=predict)
    data$predRec <- rbind(data$predRec, toAdd)

    data$selCriterion <- list(popID = popID, criterion = "pred")
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
