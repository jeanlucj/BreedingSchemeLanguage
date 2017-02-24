#'Genomic prediction
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param popID population ID to be predicted (default: the latest population)
#'@param trainingPopID population ID to be used for training a prediction model (default: all populations with phenotype data). NOTE: if sharingInfo="none" this parameter has no effect.
#'@param locations data from which locations should be used (default: all locations)
#'@param years data from which years should be used (default: all years)
#'@param sharingInfo one of "none", "markers", "pedigree".  If none, genotypic values are assumed IID. If markers or pedigree, a genomic or pedigree relationship matrix is constructed
#'
#'@return predicted values and the all information created before (list)
#'
#'@export
predictValue <- function(sEnv=simEnv, popID=NULL, trainingPopID=NULL, locations=NULL, years=NULL, sharingInfo="none"){
  parent.env(sEnv) <- environment()
  predict.func <- function(data, popID, trainingPopID, locations, years, sharingInfo){
    if (is.null(locations)) locations <- unique(data$phenoRec$loc)
    if (is.null(years)) years <- unique(data$phenoRec$year)
    phenoRec <- subset(data$phenoRec, subset=loc %in% locations & year %in% years)
    
    if (sharingInfo == "none"){
      fitIID <- lmer(formula=pValue ~ year*loc + (1|phenoGID) + (1|phenoGID:year) + (1|phenoGID:loc), data=phenoRec, weights=1/phenoRec$error)
      predict <- ranef(fitIID)$phenoGID
      predGID <- as.numeric(rownames(predict))
    }
    
    if (sharingInfo == "markers"){
      trainCandidates <- data$genoRec$hasGeno & data$genoRec$GID %in% phenoRec$phenoGID
      if(is.null(trainingPopID)){
        GID.train <- data$genoRec$GID[trainCandidates]
      }else{
        tf <- data$genoRec$popID %in% trainingPopID
        GID.train <- data$genoRec$GID[tf & trainCandidates]
      }
      mrkPos <- data$mapData$markerPos
      M <- (data$geno[GID.train*2 - 1, mrkPos] + data$geno[GID.train*2, mrkPos]) / 2
      
      # Figure out who to predict
      if (is.null(popID)) popID <- max(data$genoRec$popID)
      tf <- data$genoRec$popID %in% popID
      GID.pred <- setdiff(data$genoRec$GID[tf], GID.train)
      M <- rbind(M, (data$geno[GID.pred*2 - 1, mrkPos] + data$geno[GID.pred*2, mrkPos]) / 2)
      predGID <- c(GID.train, GID.pred)
      mt1ObsPerGID <- sum(phenoRec$phenoGID %in% predGID) > length(predGID)
      
      K <- A.mat(M)
      rownames(K) <- colnames(K) <- predGID
      kbDat <- subset(phenoRec, phenoRec$phenoGID %in% predGID)
      kbo <- kin.blup(kbDat, geno="GID", pheno="pheno", fixed=c("loc", "year"), K=K, reduce=mt1ObsPerGID, R=kbDat$error)
      predict <- kbo$g
    }
    
    if (sharingInfo == "pedigree"){
      data$aMat <- calcAmatrix(data$genoRec[, 1:3], "Outbred", data$aMat)
      
      trainCandidates <- data$genoRec$GID %in% phenoRec$phenoGID
      if(is.null(trainingPopID)){
        GID.train <- data$genoRec$GID[trainCandidates]
      }else{
        tf <- data$genoRec$popID %in% trainingPopID
        GID.train <- data$genoRec$GID[tf & trainCandidates]
      }
      if (is.null(popID)) popID <- max(data$genoRec$popID)
      tf <- data$genoRec$popID %in% popID
      GID.pred <- setdiff(data$genoRec$GID[tf], GID.train)
      predGID <- c(GID.train, GID.pred)
      mt1ObsPerGID <- sum(phenoRec$phenoGID %in% predGID) > length(predGID)
      
      K <- data$aMat[predGID, predGID]
      rownames(K) <- colnames(K) <- predGID
      kbDat <- subset(phenoRec, phenoRec$phenoGID %in% predGID)
      kbo <- kin.blup(kbDat, geno="GID", pheno="pheno", fixed=c("loc", "year"), K=K, reduce=mt1ObsPerGID, R=kbDat$error)
      predict <- kbo$g
    }
        
    if(is.null(data$predRec)){
      predNo <- 1
    } else{
      predNo <- max(data$predRec$predNo) + 1
    }
    toAdd <- data.frame(predGID=predGID, predNo=predNo, predict=predict)
    data$predRec <- rbind(data$predRec, toAdd)

    data$selCriterion <- list(popID=popID, criterion="pred")
    if (exists("totalCost", data)) data$totalCost <- data$totalCost + data$costs$predCost
    return(data)
  }#END predict.func
  with(sEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, predict.func, popID=popID, trainingPopID=trainingPopID, locations=locations, years=years, sharingInfo=sharingInfo)
      sfStop()
    }else{
      sims <- lapply(sims, predict.func, popID=popID, trainingPopID=trainingPopID, locations=locations, years=years, sharingInfo=sharingInfo)
    }
  })
}
