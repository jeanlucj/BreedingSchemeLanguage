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
      muSim$muSim <- muSim$muSim[, as.character(popID), drop=F]
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
    nrp <- 0
    if (nrow(muSim) > 1){
      g <- c(g, apply(muSim, 2, mean))
      group <- c(group, rep(sEnv$nSim + 1, nGenPlot))
      size <- c(size, rep(2, nGenPlot))
      nrp <- 1
    }
    plotData <- data.frame(g=g, popID=rep(0:(nGenPlot - 1), nrow(muSim) + nrp), size=size, col=loc, group=group, scheme=1)
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
  
  totCost <- sEnv$sims[[1]]$totalCost
  if (add){
    prevData <- try(suppressWarnings(readRDS(file=paste(addDataFileName, ".rds", sep=""))), silent=T)
    if (class(prevData) != "try-error"){
    totCost <- c(prevData$totCost, totCost)
    prevData <- prevData$plotData
    plotData$scheme <- plotData$scheme + max(prevData$scheme)
    plotData$group <- plotData$group + max(prevData$group)
    plotData <- rbind(plotData, prevData)
    }
  }
  saveRDS(list(plotData=plotData, totCost=totCost), file=paste(addDataFileName, ".rds", sep=""))
  
  mapping <- aes(x=popID, y=g, group=group)
  if (length(unique(plotData$col)) > 1) mapping <- modifyList(mapping, aes(colour=factor(col)))
  if (length(unique(plotData$size)) > 1) mapping <- modifyList(mapping, aes(size=factor(size)))
  if (length(unique(plotData$scheme)) > 1) mapping <- modifyList(mapping, aes(linetype=factor(scheme)))
  p <- ggplot(data=plotData, mapping)
  p <- p + geom_line()
  if (is.null(ymax)) {
    p <- p + ylim(min(plotData$g), max(plotData$g))
  }
  else {
    p <- p + ylim(min(plotData$g), ymax)
  }
  p <- p + labs(title="", x="Generation", y="Genetic improvement")
  
  if (length(unique(plotData$col)) > 1){
    cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    p <- p + scale_colour_manual(values=cbPalette)
    p <- p + guides(col=guide_legend("Locs"))
  } 
  if (length(unique(plotData$size)) > 1){
    p <- p + scale_size_manual(name="", values=c(0.3, 2), labels=c("Repl", "Mean"))
    p <- p + guides(size=guide_legend("Lines"))
  } 
  if (length(unique(plotData$scheme)) > 1){
    p <- p + guides(linetype=guide_legend("Scheme"))
  }
  if (!is.null(totCost)){
    p <- p + ggtitle(paste("Cost of scheme", ifelse(length(totCost) > 1, "s", ""), ": ", paste(round(totCost), collapse=", "), sep=""))
  }
  print(p)
}
