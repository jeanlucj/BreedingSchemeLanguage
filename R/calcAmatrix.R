#' Calculate an additive relationship matrix
#'
#' \code{pedigreeToAmatrix} returns an additive relationship matrix from a
#'  pedigree specified in three columns. The first column has to be the row
#'  number and sire and dam columns refer directly to rows of the pedigree.
#'
#' \code{pedigreeToAmatrix} has some functionality useful for plants.  It can
#'  handle a pedigree of doubled haploid lines. Individuals can be
#'  self-pollinated for an arbitrary number of generations.
#'
#' @param pedColumns A data.frame with four columns. The first column
#'  has to be the row number and sire and dam columns refer directly to rows
#'  of the pedigree. Parents of founders need to be set to 0 (ZERO). The row
#'  of a child has to be after (i.e. a higher row number) that of its parents.
#'  If an individual has one known and one unknown parent, set the unknown
#'  parent to 0. The fourth column indicates: If negative, the individual is
#'  a DH: the F1 is created, generates a gamete, which is then doubled.
#'  If positive, the number of generations an individual was self-pollinated
#'  after it's F1 ancestor was created (can be 0 if the individual is the F1).
#'  
#' @param aMatIn A square matrix that contains the additive relationship
#'  matrix between individuals at the beginning of the pedigree. If given,
#'  the function saves time by starting calculations after those individuals
#'  This aMatIn functionality is NOT compatible with calculating A inverse
#'
#' @return A matrix, \code{aMat}, the additive relationship matrix 
#'
calcAmatrix <- function(pedColumns, aMatIn=NULL){
  calcAmatRow <- function(pedRec){ # Function to process one row of pedColumns
    prog <- pedRec[1]
    sire <- max(pedRec[2], 0) # Founder population has negative parents
    dam <- max(pedRec[3], 0)
    progM1 <- prog - 1
    if (sire){
      sireRow <- aMat[sire, 1:progM1]
    } else{
      sireRow <- integer(progM1)
    }
    if (dam){
      damRow <- aMat[dam, 1:progM1]
    } else{
      damRow <- integer(progM1)
    }
    aMat[prog, 1:progM1] <<- (sireRow + damRow) / 2
    aMat[1:progM1, prog] <<- (sireRow + damRow) / 2
    if (pedRec[4] < 0){
      aSelf <- 2
    } else{
      if (sire > 0 & dam > 0){
        aSelf <- 1 + aMat[sire, dam] / 2
      } else{
        aSelf <- 1
      }
      if (pedRec[4] > 0){ # Number generations individual was selfed
        for (i in 1:pedRec[4]) aSelf <- 1 + aSelf / 2
      }
    }
    aMat[prog, prog] <<- aSelf
  }#END calcAmatRow

  # calculate A here
  nInd <- nrow(pedColumns)
  aMat <- matrix(0, nInd, nInd)
  if (is.null(aMatIn)){
    # the very first individual in the pedigree has to be a founder
    initSelf <- 1
    if (pedColumns[1,4] < 0) initSelf <- 2
    if (pedColumns[1,4] > 0) initSelf <- 2 - 2^(-pedColumns[1,4])
    aMatIn <- as.matrix(initSelf)
  }
  nIn <- nrow(aMatIn)
  aMat[1:nIn, 1:nIn] <- aMatIn
  start <- nIn + 1
  apply(pedColumns[start:nInd,], 1, calcAmatRow)
  rownames(aMat) <- colnames(aMat) <- pedColumns[,1]
  return(aMat)
}
