#'Create a founder population
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param nInd population size
#'
#'@return initial population informationand the all information created before (list)
#'
#'@export
initializePopulation <- function(sEnv=simEnv, nInd=100){
  parent.env(sEnv) <- environment()
  initializePopulation.func <- function(data, nInd){
    seed <- round(runif(1, 0, 1e9))
    md <- data$mapData
    
    geno <- data$founderHaps * 2 - 1
    data$founderHaps <- NULL
    geno <- geno[sample(nrow(geno), nInd*2, replace=T),]
    geno <- randomMate(popSize=nInd, geno=geno, pos=md$map$Pos)
    pedigree <- -geno$pedigree # For founders, parents will be negative
    geno <- geno$progenies
    
    M <- (geno[1:nInd*2 - 1, md$effectivePos] + geno[1:nInd*2, md$effectivePos]) / 2
    M <- scale(M, center=T, scale=F)
    mrkCenter <- attr(M, "scaled:center")
    mrkConst <- mrkCenter / 2 + 0.5; mrkConst <- 2 * crossprod(mrkConst, 1 - mrkConst)
    qtlRelMat <- tcrossprod(scale(M, center=T, scale=F)) / c(mrkConst)
    # Genetic effects. This works even if locCov is scalar
    gValue <- calcGenotypicValue(geno=geno, mapData=md)
    coef <- solve(chol(var(gValue))) %*% chol(data$varParms$locCov)
    md$effects <- md$effects %*% coef
    gValue <- gValue %*% coef
    # Year and location effects: create matrices with zero columns until phenotyped
    locEffects <- matrix(0, nrow=nInd, ncol=0)
    locEffectsI <- matrix(0, nrow=nInd, ncol=0)
    yearEffects <- matrix(0, nrow=nInd, ncol=0)
    yearEffectsI <- matrix(0, nrow=nInd, ncol=0)
    
    genoRec <- data.frame(GID=1:nInd, pedigree=pedigree, popID=0, basePopID=0, hasGeno=FALSE)
    data$mapData <- md
    data$geno <- geno; data$genoRec <- genoRec; data$gValue <- gValue
    data$locEffects <- locEffects; data$locEffectsI <- locEffectsI
    data$yearEffects <- yearEffects; data$yearEffectsI <- yearEffectsI
    data$mrkCenter <- mrkCenter; data$mrkConst <- mrkConst; data$qtlRelMat <- qtlRelMat
    return(data)
  }
  with(sEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, initializePopulation.func, nInd=nInd)
      sfStop()
    }else{
      sims <- lapply(sims, initializePopulation.func, nInd=nInd)
    }
  })
}
