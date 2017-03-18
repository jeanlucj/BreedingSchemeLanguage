#'Select individuals
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param nSelect the number of selected individuals
#'@param popID population ID to be selected (default: When random=T, the last population. When random=F, it is the last evaluated population)
#'@param random assuming random selection or selection according to their features (T: random selection, F: selection of good individuals)
#'@param type "WithinFamily", "Mass" (the default), or "OptimContrib". If Mass, all individuals are ranked against each other and the highest nSelect are taken.  If WithinFamily, individuals are ranked within half-sib (if population was randomly mated) or full-sib (if population from selfFertilize or doubledHaploid) the highest nSelect within families are taken. If "OptimContrib" all individuals are in the same pool but Meuwissen's optimum contribution is calculated using the marker-base kinship matrix.  NOTE: new calls using "OptimContrib" write over previous optimal contribution estimates.
#'
#'@return information of the selected individuals and the all information created before (list)
#'
#'@export
select <- function(sEnv=simEnv, nSelect=40, popID=NULL, random=F, type="Mass"){
  parent.env(sEnv) <- environment()
  select.func <- function(bsl, nSelect, popID, random, type){
    criterion <- bsl$selCriterion$criterion
    if(is.null(popID)){
      popID <- bsl$selCriterion$popID
      if(is.null(popID)) popID <- 0
    }
    tf <- bsl$genoRec$popID %in% popID
    GIDcan <- bsl$genoRec$GID[tf]
    if (random){
      selectedGID <- sample(GIDcan, nSelect)
    } else{
      if(substr(criterion, 1, 5) == "pheno"){
        usePheno <- subset(bsl$phenoRec, phenoGID %in% GIDcan)
        GIDcan <- intersect(GIDcan, usePheno$phenoGID)
        candValue <- by(usePheno, as.factor(usePheno$phenoGID), function(gidRec) weighted.mean(x=gidRec$pValue, w=1/gidRec$error))
      }
      if(substr(criterion, 1, 4) == "pred"){
        usePred <- subset(bsl$predRec, predGID %in% GIDcan & predNo == max(predNo))
        GIDcan <- intersect(GIDcan, usePred$predGID)
        candValue <- usePred$predict[order(usePred$predGID)]
      }
      if (type == "WithinFamily"){
        canVal <- data.frame(GID=GIDcan, val=as.numeric(candValue))
        raggedSel <- function(canVal){
          nSel <- min(nrow(canVal), nSelect)
          return(canVal$GID[order(canVal$val, decreasing=T)[1:nSel]])
        }
        selectedGID <- unlist(by(canVal, as.factor(bsl$genoRec$pedigree.1[GIDcan]), raggedSel))
      } else{
        if (type == "OptimContrib"){
          optCont <- useMeuwissenOptCont(bsl, GIDcan, candValue, nSelect)
          selectedGID <- GIDcan[optCont$optContrib > 0]
          optimContrib <- optCont$optContrib[optCont$optContrib > 0]
          bsl$selCriterion$optimContrib <- optimContrib
          bsl$selCriterion$criterion <- "OptimContrib"
          if(is.null(bsl$optContRec)){
            optContNo <- 1
          } else{
            optContNo <- max(bsl$optContRec$optContNo) + 1
          }
          toAdd <- data.frame(selectedGID, optContNo, optimContrib)
          colnames(toAdd) <- c("selectedGID", "optContNo", "optimContrib")
          bsl$optContRec <- rbind(bsl$optContRec, toAdd)
        } else{
          selectedGID <- GIDcan[order(candValue, decreasing=T)[1:nSelect]]
        }
      }
    }#END not random selection
    popID.new <- max(bsl$genoRec$popID) + 1
    bsl$genoRec$popID[bsl$genoRec$GID %in% selectedGID] <- popID.new
    if (exists("totalCost", bsl)) bsl$totalCost <- bsl$totalCost + bsl$costs$selectCost
    return(bsl)
  } #END select.func
  with(sEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, select.func, nSelect=nSelect, popID=popID, random=random, type=type)
      sfStop()
    }else{
      sims <- lapply(sims, select.func, nSelect=nSelect, popID=popID, random=random, type=type)
    }
  })
}


#' Meuwissen's optimal contribution algorithm
#'
#'@param target The target inbreeding level
#'@param varGen The kinship matrix among individuals
#'@param gebvs The selection criterion (high value is desired)
#'
#'@return vector showing the optimal contributions
#'
meuwissenOptCont <- function(target, varGen, gebvs){
  nSelCan <- length(gebvs)
  idxToOrigin <- 1:nSelCan
  anotherIter <- TRUE
  while(anotherIter){
    varGenInv <- solve(varGen)
    sumGInv <- sum(varGenInv)
    gTGInv <- t(gebvs) %*% varGenInv
    lambda02 <- (gTGInv %*% gebvs * sumGInv - sum(gTGInv)^2) / (8 * target * sumGInv - 4)
    if (lambda02 < 0){
      cat(cycle, "lambda0^2 negative:", (gTGInv %*% gebvs * sumGInv - sum(gTGInv)^2), (8 * target * sumGInv - 4), lambda02, "\n")
      lambda02 <- 1e-6
    }
    lambda0 <- c(sqrt(lambda02))
    lambda <- c(colSums(varGenInv) %*% gebvs - 2*lambda0)/sumGInv
    meuwOptCont <- varGenInv %*% (gebvs - lambda) / (2*lambda0)
    if (anotherIter <- any(meuwOptCont < 0)){
      negCont <- which.min(meuwOptCont)
      varGen <- varGen[-negCont, -negCont]
      gebvs <- gebvs[-negCont]
      idxToOrigin <- idxToOrigin[-negCont]
    }
  }
  finalMeuwOptCont <- numeric(nSelCan)
  finalMeuwOptCont[idxToOrigin] <- meuwOptCont / sum(meuwOptCont)
  return(finalMeuwOptCont)
}#END meuwissenOptCont

#' Meuwissen's optimal contribution algorithm
#'
#'@param bsl the BSL list with all objects for a simulation
#'@param gebvs The selection criterion (high value is desired) should be nInd long
#'@param selectSet should be the individuals selected by truncation selection (or whatever the "standard" selection procedure is)
#'
#'@return list with: optContrib vector showing the optimal contributions, target showing the expected inbreeding, and stdSelInb showing inbreeding under the standard method
#'
useMeuwissenOptCont <- function(bsl, GIDcan, gebvs, nSelect){
  mrkPos <- bsl$mapData$markerPos
  M <- (bsl$geno[GIDcan*2 - 1, mrkPos] + bsl$geno[GIDcan*2, mrkPos]) / 2
  nInd <- length(GIDcan)
  varGen <- A.mat(M)
  rownames(M) <- colnames(M) <- GIDcan
  selectSet <- GIDcan[order(gebvs, decreasing=T)[1:nSelect]]
  # Determine the target level of inbreeding
  cTrunc <- numeric(nInd)
  cTrunc[GIDcan %in% selectSet] <- 1 / length(selectSet) # Assumes equal contribution among selected parents under standard selection
  stdSelInb <- c(t(cTrunc) %*% varGen %*% cTrunc) / 2
  target <- stdSelInb / 2 # Hard coded to seek half the rate of inbreeding
  optCont <- try(meuwissenOptCont(target, varGen, gebvs), silent=TRUE)
  if (class(optCont) == "try-error"){
    cat(cycle, "Opt cont FAILED. Using trunc sel.", "\n")
    optCont <- cTrunc
  }
  cat(cycle, "Trunc sel inb:", stdSelInb, "Opt cont inb:", target, "Mean:", mean(gebvs), "Trunc sel gain:", c(t(cTrunc) %*% gebvs), "Opt cont gain:", c(t(optCont) %*% gebvs), "\n")
  return(list(optContrib=optCont, target=target, stdSelInb=stdSelInb))
}
