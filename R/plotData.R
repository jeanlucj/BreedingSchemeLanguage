#'Plot the results
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param ymax the maximum value of y-axis (default: the maximun value in the data)
#'@param add if T a new result will be added to an existing plot, if F a new plot will be drawn (default)
#'@param addDataFileName the name to save the summarized data for the next simulation with double-quotation, like "plot1_1". (default: "plotData")
#'@param popID vector of the population IDs you want plotted
#'
#'@return A ggplot object of the simulation results
#'
#'@export
plotData <- function(sEnv=simEnv, ymax=NULL, add=F, addDataFileName="plotData", popID=NULL){
  plotBase <- is.null(popID)
  if (plotBase) popID <- sort(unique(sEnv$sims[[1]]$genoRec$basePopID))
  nLoc <- ncol(sEnv$sims[[1]]$gValue)
  
  getMeans <- function(bsl, loc){
    if (plotBase) pID <- bsl$genoRec$basePopID
    else pID <- bsl$genoRec$popID
    return(tapply(bsl$gValue[,loc], pID, mean))
  }
  muSimByLoc <- lapply(1:nLoc, function(loc) list(muSim=t(sapply(sEnv$sims, getMeans, loc=loc)), loc=loc))
  
  if (class(popID) == "list"){
    pID <- sEnv$sims[[1]]$genoRec$popID
    popSizes <- tapply(pID, pID, length)
    modifyMSBL <- function(muSim){
      ms <- muSim$muSim
      muSim$muSim <- sapply(popID, function(popVec) apply(ms, 1, function(vec) weighted.mean(vec[as.character(popVec)], popSizes[as.character(popVec)])))
      return(muSim)
    }
  } else{
    modifyMSBL <- function(muSim){
      muSim$muSim <- muSim$muSim[, as.character(popID)]
      return(muSim)
    }
  }
  muSimByLoc <- lapply(muSimByLoc, modifyMSBL)
  
  makeDF <- function(muSim){
    loc <- muSim$loc
    muSim <- muSim$muSim
    muSim <- muSim - muSim[, 1]
    g <- NULL
    group <- NULL
    size <- NULL
    nGenPlot <- length(popID)
    for(sim in 1:nrow(muSim)){
      g <- c(g, muSim[sim, ])
      group <- c(group, rep(sim, nGenPlot))
      size <- c(size, rep(1, nGenPlot))
    }
    g <- c(g, apply(muSim, 2, mean))
    group <- c(group, rep(sEnv$nSim + 1, nGenPlot))
    size <- c(size, rep(2, nGenPlot))
    plotData <- data.frame(g=g, popID=rep(0:(nGenPlot - 1), nrow(muSim) + 1), size=size, col=loc, group=group, scheme=1)
  }#END makeDF
  muSimByLoc <- lapply(muSimByLoc, makeDF)
  
  plotData <- NULL
  maxGroup <- 0
  for (loc in 1:nLoc){
    toAdd <- muSimByLoc[[loc]]
    toAdd$group <- toAdd$group + maxGroup
    maxGroup <- max(toAdd$group)
    plotData <- rbind(plotData, toAdd)
  }
  
  if (add){
    prevData <- readRDS(file=paste(addDataFileName, ".rds", sep=""))
    plotData$scheme <- plotData$scheme + max(prevData$scheme)
    plotData$group <- plotData$group + max(prevData$group)
    plotData <- rbind(plotData, prevData)
  }
  saveRDS(plotData, file=paste(addDataFileName, ".rds", sep=""))
  p <- ggplot(data=plotData, aes(x=popID, y=g))
  p <- p + geom_line(aes(size=factor(size), colour=factor(col), linetype=factor(scheme), group=factor(group)))
  if (is.null(ymax)) {
    p <- p + ylim(min(plotData$g), max(plotData$g))
  }
  else {
    p <- p + ylim(min(plotData$g), ymax)
  }
  # The palette with black:
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # To use for line and point colors, add
  p <- p + scale_colour_manual(values=cbPalette)
  p <- p + scale_size_manual(name="", values=c(0.3, 2), labels=c("Repl", "Mean"))
  p <- p + labs(title="", x="Generation", y="Genetic improvement")
  p <- p + guides(col=guide_legend("Locs"))
  p <- p + guides(size=guide_legend("Lines"))
  p <- p + guides(linetype=guide_legend("Scheme"))
  if (exists("totalCost", sEnv$sims[[1]])){
    p <- p + ggtitle(paste("Total scheme cost", round(sEnv$sims[[1]]$totalCost)))
  }
  print(p)
}
