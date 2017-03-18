#'Select individuals
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param nSelect the number of selected individuals
#'@param popID population ID to be selected (default: When random=T, the last population. When random=F, it is the last evaluated population)
#'@param random assuming random selection or selection according to their features (T: random selection, F: selection of good individuals)
#'@param type "WithinFamily" or "Mass" (default: Mass). If Mass, all individuals are ranked against each other and the highest nSelect are taken.  If WithinFamily, individuals are ranked within half-sib (if population was randomly mated) or full-sib (if population from selfFertilize or doubledHaploid) the highest nSelect within families are taken.
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
        GIDcan <- intersect(GIDcan, bsl$phenoRec$phenoGID)
        usePheno <- subset(bsl$phenoRec, phenoGID %in% GIDcan)
        candValue <- by(usePheno, as.factor(usePheno$phenoGID), function(gidRec) weighted.mean(x=gidRec$pValue, w=1/gidRec$error))
      }
      if(substr(criterion, 1, 4) == "pred"){
        GIDcan <- intersect(GIDcan, bsl$predRec$predGID)
        usePred <- bsl$predRec[bsl$predRec$predGID %in% GIDcan & bsl$predRec$predNo == max(bsl$predRec$predNo),]
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
        selectedGID <- GIDcan[order(candValue, decreasing=T)[1:nSelect]]
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
