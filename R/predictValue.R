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
  predictValue.func <- function(bsl, popID, trainingPopID, locations, years, sharingInfo){
    if (is.null(locations)) locations <- unique(bsl$phenoRec$loc)
    if (is.null(years)) years <- unique(bsl$phenoRec$year)
    phenoRec <- subset(bsl$phenoRec, subset=loc %in% locations & year %in% years)
    ########################################################
    # Consider all individuals to be IID
    if (sharingInfo == "none"){
      mt1ObsPerGID <- length(unique(phenoRec$phenoGID)) < nrow(phenoRec)
      if (mt1ObsPerGID){ # More than one observation per GID: run a model
        fmla <- "pValue ~ (1|phenoGID)"
        if (length(unique(phenoRec$year))) fmla <- paste(fmla, "+ year + (1|phenoGID:year)")
        if (length(unique(phenoRec$loc))) fmla <- paste(fmla, "+ loc + (1|phenoGID:loc)")
        phenoRec$phenoGID <- factor(phenoRec$phenoGID)
        phenoRec$loc <- factor(phenoRec$loc)
        phenoRec$year <- factor(phenoRec$year)
        fitIID <- lmer(formula=as.formula(fmla), data=phenoRec, weights=1/phenoRec$error)
        predict <- ranef(fitIID)$phenoGID
        predGID <- as.numeric(rownames(predict))
      } else{ # Only one observation per GID
        predict <- phenoRec$pValue
        predGID <- as.numeric(phenoRec$phenoGID)
      }
    }
    ########################################################
    # Use markers to determine individual relatedness
    if (sharingInfo == "markers"){
      trainCandidates <- bsl$genoRec$hasGeno & bsl$genoRec$GID %in% phenoRec$phenoGID
      if(is.null(trainingPopID)){
        GID.train <- bsl$genoRec$GID[trainCandidates]
      }else{
        tf <- bsl$genoRec$popID %in% trainingPopID
        GID.train <- bsl$genoRec$GID[tf & trainCandidates]
      }
      mrkPos <- bsl$mapData$markerPos
      M <- (bsl$geno[GID.train*2 - 1, mrkPos] + bsl$geno[GID.train*2, mrkPos]) / 2
      
      # Figure out who to predict
      if (is.null(popID)) popID <- max(bsl$genoRec$popID)
      tf <- bsl$genoRec$popID %in% popID
      GID.pred <- setdiff(bsl$genoRec$GID[tf], GID.train)
      M <- rbind(M, (bsl$geno[GID.pred*2 - 1, mrkPos] + bsl$geno[GID.pred*2, mrkPos]) / 2)
      predGID <- c(GID.train, GID.pred)
      mt1ObsPerGID <- sum(phenoRec$phenoGID %in% predGID) > length(predGID)
      
      K <- A.mat(M)
      rownames(K) <- colnames(K) <- predGID
      kbDat <- subset(phenoRec, phenoRec$phenoGID %in% predGID)
      kbo <- kin.blup(kbDat, geno="phenoGID", pheno="pValue", fixed=c("loc", "year"), K=K, reduce=mt1ObsPerGID, R=kbDat$error)
      predict <- kbo$g
    }
    ########################################################
    # Use pedigree to determine individual relatedness
    if (sharingInfo == "pedigree"){
      bsl$aMat <- calcAmatrix(bsl$genoRec[, 1:4], bsl$aMat)
      
      trainCandidates <- bsl$genoRec$GID %in% phenoRec$phenoGID
      if(is.null(trainingPopID)){
        GID.train <- bsl$genoRec$GID[trainCandidates]
      }else{
        tf <- bsl$genoRec$popID %in% trainingPopID
        GID.train <- bsl$genoRec$GID[tf & trainCandidates]
      }
      if (is.null(popID)) popID <- max(bsl$genoRec$popID)
      tf <- bsl$genoRec$popID %in% popID
      GID.pred <- setdiff(bsl$genoRec$GID[tf], GID.train)
      predGID <- c(GID.train, GID.pred)
      mt1ObsPerGID <- sum(phenoRec$phenoGID %in% predGID) > length(predGID)
      
      K <- bsl$aMat[predGID, predGID]
      rownames(K) <- colnames(K) <- predGID
      kbDat <- subset(phenoRec, phenoRec$phenoGID %in% predGID)
      kbo <- kin.blup(kbDat, geno="phenoGID", pheno="pValue", fixed=c("loc", "year"), K=K, reduce=mt1ObsPerGID, R=kbDat$error)
      predict <- kbo$g
    }
        
    if(is.null(bsl$predRec)){
      predNo <- 1
    } else{
      predNo <- max(bsl$predRec$predNo) + 1
    }
    toAdd <- data.frame(predGID, predNo, predict)
    colnames(toAdd) <- c("predGID", "predNo", "predict")
    bsl$predRec <- rbind(bsl$predRec, toAdd)

    bsl$selCriterion <- list(popID=popID, criterion="pred")
    if (exists("totalCost", bsl)) bsl$totalCost <- bsl$totalCost + bsl$costs$predCost
    return(bsl)
  }#END predict.func
  with(sEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, predictValue.func, popID=popID, trainingPopID=trainingPopID, locations=locations, years=years, sharingInfo=sharingInfo)
      sfStop()
    }else{
      sims <- lapply(sims, predictValue.func, popID=popID, trainingPopID=trainingPopID, locations=locations, years=years, sharingInfo=sharingInfo)
    }
  })
}
