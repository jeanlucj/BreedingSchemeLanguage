#'Select individuals
#'
#'@param simEnv an environment that BSL statements operate on
#'@param nSelect the number of selected individuals
#'@param popID population ID to be selected (default: When random=T, the last population. When random=F, it is the last evaluated population)
#'@param random assuming random selection or selection according to their features (T: random selection, F: selection of good individuals)
#'
#'@return information of the selected individuals and the all information created before (list)
#'
#'@export
select <- function(simEnv, nSelect=40, popID=NULL, random=F){
  parent.env(simEnv) <- environment()
  select.func <- function(data, nSelect, popID, random=FALSE){
    criterion <- data$selCriterion$criterion
    if(is.null(popID)){
      popID <- data$selCriterion$popID
      if(is.null(popID)) popID <- 0
    }
    tf <- data$genoRec$popID %in% popID
    GIDcan <- data$genoRec$GID[tf]
    if (random){
      selectedGID <- sample(GIDcan, nSelect)
    } else{
      if(substr(criterion, 1, 5) == "pheno"){
        usePheno <- data$phenoRec[phenoRec$phenoGID %in% GIDcan,]
        candValue <- by(usePheno, as.factor(usePheno$phenoGID), function(gidRec) weighted.mean(x=gidRec$pValue, w=1/gidRec$error))
      }else{
        if(substr(criterion, 1, 4) == "pred"){
          usePred <- data$predRec[data$predRec$predGID %in% GIDcan & data$predRec$predNo == max(data$predRec$predNo),]
          candValue <- usePred$predict[order(usePred$predGID)]
        }else{
          stop("Please define selection criterion in correct way!")
        }
      }
      selectedGID <- GIDcan[order(candValue, decreasing=T)[1:nSelect]]
    }#END not random selection
    popID.new <- max(data$genoRec$popID) + 1
    data$genoRec$popID[data$genoRec$GID %in% selectedGID] <- popID.new
    return(data)
  } #END select.func
  with(simEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, select.func, nSelect=nSelect, popID=popID, random=random)
      sfStop()
    }else{
      sims <- lapply(sims, select.func, nSelect=nSelect, popID=popID, random=random)
    }
  })
}
