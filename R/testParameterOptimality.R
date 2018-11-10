#' Function to return the optimality of a parameter vector for a breeding scheme given a simulation environment
#'
#'@param sEnv the environment that BSL functions operate in. This environment should be fresh from defineSpecies, defineVariances, and defineCosts. If NULL, the default \code{simEnv} is attempted
#'@param schemeFileName source file that holds the script of the breeding scheme. The scheme should return a value to be maximized in a variable called \code{objective}. If you want to use a non-default simulation environment, you need to plan for that in the schemeFile by giving sEnv as the first parameter in all of the BSL functions
#'@param parmList (preferably) named list with values of parameters characterizing the breeding scheme
#'@param objectiveFunc a function that can be applied to a BSL simulation environment to return a value showing the realized performance of the breeding scheme
#'@param budget the maximum budget that is allowed for the breeding scheme
#'
#'@return Two-object list: \code{objective} is the objective function value of the breeding scheme given the simulation environment and the parameter vector and \code{totalCost} is the budget used by the breeding scheme. If the parameter vector causes the scheme to exceed the given budget, \code{objective}=NA
#'
#'@export
testParameterOptimality <- function(sEnv=NULL, schemeFileName, parmList, objectiveFunc, budget=1000, runChecks=T){
  if(is.null(sEnv)){
    if(exists("simEnv", .GlobalEnv)){
      sEnv <- get("simEnv", .GlobalEnv)
    } else{
      stop("No simulation environment was passed")
    }
  } 

  sEnv$onlyCost <- TRUE
  source(schemeFileName, local=TRUE)
  if (runChecks){
    # Check if there is g x year but only one year was run
    if (sEnv$costs$nYear == 1 & sEnv$sims[[1]]$varParms$gByYearVar > 0){
      print("Genotype by year interaction greater than zero but only one year was simulated. This situation causes divergence between performance (G + GxY) and what matters (G).")
    }
    # Check if there is g x loc but only one loc was run
    if (sEnv$costs$nLoc == 1 & sEnv$sims[[1]]$varParms$gByLocVar > 0){
      print("Genotype by location interaction greater than zero but only one location was simulated. This situation causes divergence between performance (G + GxL) and what matters (G).")
    }
  }
  totalCost <- sEnv$totalCost
  if (totalCost > budget) return(list(objective=NA, totalCost=totalCost))
  sEnv$onlyCost <- FALSE
  source(schemeFileName, local=TRUE)
  return(list(objective=objectiveFunc(sEnv), totalCost=totalCost))
}
