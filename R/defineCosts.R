#' Define the costs that go into breeding
#' Default for some costs is zero because they probably belong to fixed costs
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param phenoCost named vector: names are the plot types and contents are the cost for one plot of that type (default: Standard=1; if plotTypes have been defined in defineVariances and phenoCost is not specified, all plotTypes are given a cost of 1)
#'@param genoCost scalar: cost to genotype one individual default (0.25)
#'@param crossCost scalar: cost of creating a new individual from a cross (1)
#'@param selfCost scalar: the cost of creating a selfed seed (1)
#'@param doubHapCost scalar: the cost of creating a doubled haploid seed (5)
#'@param predCost scalar: the cost of running the analysis to make predictions (0)
#'@param selectCost scalar: the cost of running the analysis to do selection (0)
#'@param locCost scalar: the cost of maintaining a location for a year (0)
#'@param yearCost scalar: the cost of program upkeep for a year (0)
#'
#'@return modifies the list sims in environment sEnv by adding cost parameters allowing the total cost of the simulated scheme to be calculated
#'
#'@export
defineCosts <- function(sEnv=NULL, phenoCost=NULL, genoCost=0.25, crossCost=1, selfCost=1, doubHapCost=5, predCost=0, selectCost=0, locCost=0, yearCost=0){
  
  if(is.null(sEnv)){
    if(exists("simEnv", .GlobalEnv)){
      sEnv <- get("simEnv", .GlobalEnv)
    } else{
      stop("No simulation environment was passed")
    }
  } 
  
  bsl <- sEnv$sims[[1]]
  if (is.null(phenoCost)){
    phenoCost <- rep(1, length(bsl$varParms$plotTypeErrVars))
    names(phenoCost) <- names(bsl$varParms$plotTypeErrVars)
  }
  sEnv$costs <- list(phenoCost=phenoCost, genoCost=genoCost, crossCost=crossCost, selfCost=selfCost, doubHapCost=doubHapCost, predCost=predCost, selectCost=selectCost, locCost=locCost, yearCost=yearCost)
  if (bsl$varParms$randLoc) sEnv$totalCost <- 0 else sEnv$totalCost <- locCost * ncol(bsl$varParms$locCov)
}
