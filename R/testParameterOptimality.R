#' Function to return the optimality of a parameter vector for a breeding scheme given a simulation environment
#'
#'@param sEnv the environment that BSL functions operate in. This environment should be fresh from defineSpecies, defineVariances, and defineCosts. If NULL, the default \code{simEnv} is attempted
#'@param schemeFileName source file that holds the script of the breeding scheme. The scheme should return a value to be maximized in a variable called \code{objective}. If you want to use a non-default simulation environment, you need to plan for that in the schemeFile by giving sEnv as the first parameter in all of the BSL functions
#'@param parmVec vector with values of parameters characterizing the breeding scheme
#'@param budget the maximum budget that is allowed for the breeding scheme
#'
#'@return Optimality of the breeding scheme given the simulation environment and the parameter vector. If parameter vector causes the scheme to exceed the budget, NA is returned
#'
#'@export
testParameterOptimality <- function(sEnv=NULL, schemeFileName, parmVec, budget){

  # Uh oh: I am going to have to completely revamp the computation for
  # onlyCost because now it does not figure out, for example, 
  # how many individuals would be genotyped and phenotyped
  sEnv$onlyCost <- TRUE
  source(schemeFileName)
  if (sEnv$totalCost > budget) return(NA)
  sEnv$onlyCost <- FALSE
  source(schemeFileName)
  return(objective)
}
