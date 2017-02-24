#' Define the costs that go into breeding
#' Default for some costs is zero because they probably belong to fixed costs
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param phenoCost two-column lookup table with error variance in first column and cost in second column (default: 1, 1)
#'@param genoCost scalar: cost to genotype one individual default (0.25)
#'@param crossCost scalar: cost of creating a new individual from a cross (1)
#'@param selfCost scalar: the cost of creating a selfed seed (1)
#'@param doubHapCost scalar: the cost of creating a doubled haploid seed (5)
#'@param predCost scalar: the cost of running the analysis to make predictions (0)
#'@param selectCost scalar: the cost of running the analysis to do selection (0)
#'@param locCost scalar: the cost of maintaining a location for a year (0)
#'@param yearCost scalar: the cost of program upkeep for a year (0)
#'
#'@return Species information and input values for the simulation (list)
#'
#'@export
defineCosts <- function(sEnv=simEnv, phenoCost=data.frame(1,1), genoCost=0.25, crossCost=1, selfCost=1, doubHapCost=5, predCost=0, selectCost=0, locCost=0, yearCost=0){
  parent.env(sEnv) <- environment()
  cost.func <- function(data){
    colnames(phenoCost) <- c("error", "cost")
    data$costs <- list(phenoCost=phenoCost, genoCost=genoCost, crossCost=crossCost, selfCost=selfCost, doubHapCost=doubHapCost, predCost=predCost, selectCost=selectCost, locCost=locCost, yearCost=yearCost)
    data$totalCost <- 0
    return(data)
  }
  with(sEnv, {
    # This is too fast to want to parallelize
    sims <- lapply(sims, cost.func)
  })
}
  