#'Genomic prediction
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param popID population ID to be predicted (default: the latest population)
#'@param trainingPopID population ID to be used for training a prediction model (default: all populations with phenotype data). NOTE: if sharingInfo="none" this parameter has no effect.
#'@param locations data from which locations should be used (default: all locations)
#'@param years data from which years should be used (default: all years)
#'@param sharingInfo one of "none", "markers", "pedigree".  If none, genotypic values are assumed IID. If markers or pedigree, a genomic or pedigree relationship matrix is constructed
#'
#'@importFrom stats as.formula
#'
#'@return predicted values and the all information created before (list)
#'
#'@export
predictValue <- function(sEnv=simEnv, popID=NULL, trainingPopID=NULL, locations=NULL, years=NULL, sharingInfo=NULL){
  parent.env(sEnv) <- environment()
  predictValue.func <- function(bsl, popID, trainingPopID, locations, years, sharingInfo){
    if (is.null(popID)) popID <- max(bsl$genoRec$popID)
    if (is.null(sharingInfo)) sharingInfo <- bsl$selCriterion$sharing
    if (is.null(locations)) locations <- unique(bsl$phenoRec$loc)
    if (is.null(years)) years <- unique(bsl$phenoRec$year)
    phenoRec <- subset(bsl$phenoRec, subset = bsl$phenoRec$loc %in% locations & bsl$phenoRec$year %in% years)
    ########################################################

    # Use kin.blup4.4 that has reduce and R features
    ########################################
    kin.blup4.4 <- function(data,geno,pheno,GAUSS=FALSE,K=NULL,fixed=NULL,covariate=NULL,PEV=FALSE,n.core=1,theta.seq=NULL,reduce=FALSE,R=NULL)  {
        
        make.full <- function(X) {
            svd.X <- svd(X)
            r <- max(which(svd.X$d>1e-8))
            return(as.matrix(svd.X$u[,1:r]))
        }
        
        names <- colnames(data)	
        ypos <- match(pheno,names)
        if (is.na(ypos)) {
            stop("Phenotype name does not appear in data.")
        } else {
            y <- data[,ypos]
        }
        
        if (!is.null(R)&(length(R)!=length(y))) {
            stop("Length of R does not equal length of y")
        }
        
        not.miss <- which(!is.na(y))
        if (length(not.miss)<length(y)) {
            data <- data[not.miss,]
            y <- y[not.miss]
            if (!is.null(R)) {R <- R[not.miss]}
        } 
        n <- length(y)
        
        X <- matrix(1,n,1)		
        if (!is.null(fixed)) {
            p <- length(fixed)
            for (i in 1:p) {
                xpos <- match(fixed[i],names)
                xx <- factor(data[,xpos])	
                if (length(unique(xx)) > 1) {X <- cbind(X,model.matrix(~x-1,data.frame(x=xx)))}
            }
        }
        if (!is.null(covariate)) {
            p <- length(covariate)
            for (i in 1:p) {
                xpos <- match(covariate[i],names)
                X <- cbind(X,data[,xpos])
            }
        }	
        
        gid.pos <- match(geno,names)
        if (is.na(gid.pos)) {stop("Genotype name does not appear in data.")}
        
        not.miss.gid <- as.character(unique(data[,gid.pos]))
        if (is.null(K)) {
            if (reduce) {print("reduce=TRUE is not valid for independent genotypes. Proceeding without reduction.")}
            gid <- not.miss.gid
            v <- length(gid)
            Z <- matrix(0,n,v)
            colnames(Z) <- gid
            Z[cbind(1:n,match(data[,gid.pos],gid))] <- 1
            X2 <- make.full(X)
            ans <- rrBLUP::mixed.solve(y=y,X=X2,Z=Z,SE=PEV)
            resid <- y-X2%*%ans$beta-Z%*%ans$u
            if (PEV) {
                return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u,PEV=ans$u.SE^2,resid=resid))
            } else {
                return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u,resid=resid))
            }
        } else {
            
            if (class(K)=="dist") {K <- as.matrix(K)}
            gid <- rownames(K)
            ix.pheno <- match(not.miss.gid,gid)
            miss.pheno.gid <- which(is.na(ix.pheno))
            if (length(miss.pheno.gid)>0) {
                stop(paste("The following lines have phenotypes but no genotypes:",paste(not.miss.gid[miss.pheno.gid],collapse=" ")))
            }
            miss.gid <- setdiff(gid,not.miss.gid)
            ix <- c(ix.pheno,match(miss.gid,gid))
            K <- K[ix,ix]
            v <- length(not.miss.gid)
            Z <- matrix(0,n,v)
            Z[cbind(1:n,match(data[,gid.pos],not.miss.gid))] <- 1
            
            if (!is.null(R)) {
                sqrt.R <- sqrt(R)
                X2 <- X/sqrt.R
                y2 <- y/sqrt.R
                Z2 <- Z/sqrt.R
            } else {
                X2 <- X
                y2 <- y
                Z2 <- Z
            }
            
            if ((n > v)&(reduce)) {
                #transform
                w <- sqrt(diag(crossprod(Z2)))
                X2 <- make.full(crossprod(Z2,X2)/w)
                y2 <- crossprod(Z2,y2)/w
                Z2 <- cbind(diag(w),matrix(0,v,nrow(K)-v))
                reduced <- TRUE
            } else {
                X2 <- make.full(X2)
                Z2 <- cbind(Z2,matrix(0,n,nrow(K)-v))			
                reduced <- FALSE
            }
            
            rm(X,Z,y)
            
            if (!GAUSS) {
                ans <- rrBLUP::mixed.solve(y=y2,X=X2,Z=Z2,K=K,SE=PEV)
                ix <- match(gid,rownames(ans$u))
                if (reduced) {
                    resid <- NULL
                } else {			
                    resid <- y2-X2%*%ans$beta-Z2%*%ans$u
                }
                if (PEV) {
                    return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u[ix],PEV=ans$u.SE[ix]^2,resid=resid))
                } else {
                    return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u[ix],resid=resid))
                }
                
            } else {
                
                if (is.null(theta.seq)) {
                    theta <- setdiff(seq(0,max(K),length.out=11),0)
                } else {
                    theta <- theta.seq
                }
                n.profile <- length(theta)		
                ms.fun <- function(theta) {
                    soln <- list()
                    n.t <- length(theta)
                    for (i in 1:n.t) {
                        soln[[i]] <- rrBLUP::mixed.solve(y=y2,X=X2,Z=Z2,K=exp(-(K/theta[i])^2),SE=PEV)
                    }
                    return(soln)
                }
                
                if (n.core > 1) {
                    library(parallel)
                    it <- split(theta,factor(cut(theta,n.core,labels=FALSE)))
                    soln <- unlist(mclapply(it,ms.fun,mc.cores=n.core),recursive=FALSE)
                } else {
                    soln <- ms.fun(theta)
                }      
                
                LL <- rep(0,n.profile)
                for (i in 1:n.profile) {LL[i] <- soln[[i]]$LL}
                ans <- soln[[which.max(LL)]]	
                profile <- cbind(theta,LL)
                ix <- match(gid,rownames(ans$u))
                
                if (reduced) {
                    resid <- NULL
                } else {			
                    resid <- y2-X2%*%ans$beta-Z2%*%ans$u
                }
                
                if (PEV) {
                    return(list(Vg=ans$Vu,Ve=ans$Ve,profile=profile,g=ans$u[ix],PEV=ans$u.SE[ix]^2,resid=resid))
                } else {
                    return(list(Vg=ans$Vu,Ve=ans$Ve,profile=profile,g=ans$u[ix],resid=resid))
                }	
                
            } #else GAUSS
        } #else is.null(K)
    } #kin.blup
    
    ########################################################
    # Consider all individuals to be IID
    if (sharingInfo == "none"){
      # Figure out who to predict
      tf <- bsl$genoRec$popID %in% popID
      predGID <- bsl$genoRec$GID[tf]

      kbDat <- subset(phenoRec, phenoRec$phenoGID %in% predGID)
      kbo <- kin.blup4.4(kbDat, geno="phenoGID", pheno="pValue", fixed=c("loc", "year"), R=kbDat$error)
      predict <- kbo$g
      predGID <- as.integer(names(predict))
    }
    
    # Use kin.blup4.4 that has reduce and R features
    ########################################
    kin.blup4.4 <- function(data,geno,pheno,GAUSS=FALSE,K=NULL,fixed=NULL,covariate=NULL,PEV=FALSE,n.core=1,theta.seq=NULL,reduce=FALSE,R=NULL)  {
        
        make.full <- function(X) {
            svd.X <- svd(X)
            r <- max(which(svd.X$d>1e-8))
            return(as.matrix(svd.X$u[,1:r]))
        }
        
        names <- colnames(data)	
        ypos <- match(pheno,names)
        if (is.na(ypos)) {
            stop("Phenotype name does not appear in data.")
        } else {
            y <- data[,ypos]
        }
        
        if (!is.null(R)&(length(R)!=length(y))) {
            stop("Length of R does not equal length of y")
        }
        
        not.miss <- which(!is.na(y))
        if (length(not.miss)<length(y)) {
            data <- data[not.miss,]
            y <- y[not.miss]
            if (!is.null(R)) {R <- R[not.miss]}
        } 
        n <- length(y)
        
        X <- matrix(1,n,1)		
        if (!is.null(fixed)) {
            p <- length(fixed)
            for (i in 1:p) {
                xpos <- match(fixed[i],names)
                xx <- factor(data[,xpos])	
                if (length(unique(xx)) > 1) {X <- cbind(X,model.matrix(~x-1,data.frame(x=xx)))}
            }
        }
        if (!is.null(covariate)) {
            p <- length(covariate)
            for (i in 1:p) {
                xpos <- match(covariate[i],names)
                X <- cbind(X,data[,xpos])
            }
        }	
        
        gid.pos <- match(geno,names)
        if (is.na(gid.pos)) {stop("Genotype name does not appear in data.")}
        
        not.miss.gid <- as.character(unique(data[,gid.pos]))
        if (is.null(K)) {
            if (reduce) {print("reduce=TRUE is not valid for independent genotypes. Proceeding without reduction.")}
            gid <- not.miss.gid
            v <- length(gid)
            Z <- matrix(0,n,v)
            colnames(Z) <- gid
            Z[cbind(1:n,match(data[,gid.pos],gid))] <- 1
            X2 <- make.full(X)
            ans <- rrBLUP::mixed.solve(y=y,X=X2,Z=Z,SE=PEV)
            resid <- y-X2%*%ans$beta-Z%*%ans$u
            if (PEV) {
                return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u,PEV=ans$u.SE^2,resid=resid))
            } else {
                return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u,resid=resid))
            }
        } else {
            
            if (class(K)=="dist") {K <- as.matrix(K)}
            gid <- rownames(K)
            ix.pheno <- match(not.miss.gid,gid)
            miss.pheno.gid <- which(is.na(ix.pheno))
            if (length(miss.pheno.gid)>0) {
                stop(paste("The following lines have phenotypes but no genotypes:",paste(not.miss.gid[miss.pheno.gid],collapse=" ")))
            }
            miss.gid <- setdiff(gid,not.miss.gid)
            ix <- c(ix.pheno,match(miss.gid,gid))
            K <- K[ix,ix]
            v <- length(not.miss.gid)
            Z <- matrix(0,n,v)
            Z[cbind(1:n,match(data[,gid.pos],not.miss.gid))] <- 1
            
            if (!is.null(R)) {
                sqrt.R <- sqrt(R)
                X2 <- X/sqrt.R
                y2 <- y/sqrt.R
                Z2 <- Z/sqrt.R
            } else {
                X2 <- X
                y2 <- y
                Z2 <- Z
            }
            
            if ((n > v)&(reduce)) {
                #transform
                w <- sqrt(diag(crossprod(Z2)))
                X2 <- make.full(crossprod(Z2,X2)/w)
                y2 <- crossprod(Z2,y2)/w
                Z2 <- cbind(diag(w),matrix(0,v,nrow(K)-v))
                reduced <- TRUE
            } else {
                X2 <- make.full(X2)
                Z2 <- cbind(Z2,matrix(0,n,nrow(K)-v))			
                reduced <- FALSE
            }
            
            rm(X,Z,y)
            
            if (!GAUSS) {
                ans <- rrBLUP::mixed.solve(y=y2,X=X2,Z=Z2,K=K,SE=PEV)
                ix <- match(gid,rownames(ans$u))
                if (reduced) {
                    resid <- NULL
                } else {			
                    resid <- y2-X2%*%ans$beta-Z2%*%ans$u
                }
                if (PEV) {
                    return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u[ix],PEV=ans$u.SE[ix]^2,resid=resid))
                } else {
                    return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u[ix],resid=resid))
                }
                
            } else {
                
                if (is.null(theta.seq)) {
                    theta <- setdiff(seq(0,max(K),length.out=11),0)
                } else {
                    theta <- theta.seq
                }
                n.profile <- length(theta)		
                ms.fun <- function(theta) {
                    soln <- list()
                    n.t <- length(theta)
                    for (i in 1:n.t) {
                        soln[[i]] <- rrBLUP::mixed.solve(y=y2,X=X2,Z=Z2,K=exp(-(K/theta[i])^2),SE=PEV)
                    }
                    return(soln)
                }
                
                if (n.core > 1) {
                    library(parallel)
                    it <- split(theta,factor(cut(theta,n.core,labels=FALSE)))
                    soln <- unlist(mclapply(it,ms.fun,mc.cores=n.core),recursive=FALSE)
                } else {
                    soln <- ms.fun(theta)
                }      
                
                LL <- rep(0,n.profile)
                for (i in 1:n.profile) {LL[i] <- soln[[i]]$LL}
                ans <- soln[[which.max(LL)]]	
                profile <- cbind(theta,LL)
                ix <- match(gid,rownames(ans$u))
                
                if (reduced) {
                    resid <- NULL
                } else {			
                    resid <- y2-X2%*%ans$beta-Z2%*%ans$u
                }
                
                if (PEV) {
                    return(list(Vg=ans$Vu,Ve=ans$Ve,profile=profile,g=ans$u[ix],PEV=ans$u.SE[ix]^2,resid=resid))
                } else {
                    return(list(Vg=ans$Vu,Ve=ans$Ve,profile=profile,g=ans$u[ix],resid=resid))
                }	
                
            } #else GAUSS
        } #else is.null(K)
    } #kin.blup
    ##########################
    
    ########################################################
    # Use markers to determine individual relatedness
    if (sharingInfo == "markers"){
      trainCandidates <- bsl$genoRec$hasGeno & bsl$genoRec$GID %in% phenoRec$phenoGID
      if (is.null(trainingPopID)){
        GID.train <- bsl$genoRec$GID[trainCandidates]
      } else{
        tf <- bsl$genoRec$popID %in% trainingPopID
        GID.train <- bsl$genoRec$GID[tf & trainCandidates]
      }
      mrkPos <- bsl$mapData$markerPos
      M <- (bsl$geno[GID.train*2 - 1, mrkPos] + bsl$geno[GID.train*2, mrkPos]) / 2
      
      # Figure out who to predict
      tf <- bsl$genoRec$popID %in% popID
      GID.pred <- setdiff(bsl$genoRec$GID[tf], GID.train)
      M <- rbind(M, (bsl$geno[GID.pred*2 - 1, mrkPos] + bsl$geno[GID.pred*2, mrkPos]) / 2)
      predGID <- c(GID.train, GID.pred)
      mt1ObsPerGID <- sum(phenoRec$phenoGID %in% predGID) > length(predGID)
      
      K <- rrBLUP::A.mat(M)
      rownames(K) <- colnames(K) <- predGID
      kbDat <- subset(phenoRec, phenoRec$phenoGID %in% predGID)
      kbo <- kin.blup4.4(kbDat, geno="phenoGID", pheno="pValue", fixed=c("loc", "year"), K=K, reduce=mt1ObsPerGID, R=kbDat$error)
      predict <- kbo$g
      predGID <- as.integer(names(predict))
    }
    ########################################################
    # Use pedigree to determine individual relatedness
    if (sharingInfo == "pedigree"){
      bsl$aMat <- calcAmatrix(bsl$genoRec[, 1:4], bsl$aMat)
      
      trainCandidates <- bsl$genoRec$GID %in% phenoRec$phenoGID
      if(is.null(trainingPopID)){
        GID.train <- bsl$genoRec$GID[trainCandidates]
      }else{
        tf <- bsl$genoRec$popID %in% trainingPopID
        GID.train <- bsl$genoRec$GID[tf & trainCandidates]
      }
      tf <- bsl$genoRec$popID %in% popID
      GID.pred <- setdiff(bsl$genoRec$GID[tf], GID.train)
      predGID <- c(GID.train, GID.pred)
      mt1ObsPerGID <- sum(phenoRec$phenoGID %in% predGID) > length(predGID)
      
      K <- bsl$aMat[predGID, predGID]
      rownames(K) <- colnames(K) <- predGID
      kbDat <- subset(phenoRec, phenoRec$phenoGID %in% predGID)
      kbo <- kin.blup4.4(kbDat, geno="phenoGID", pheno="pValue", fixed=c("loc", "year"), K=K, reduce=mt1ObsPerGID, R=kbDat$error)
      predict <- kbo$g
      predGID <- as.integer(names(predict))
    }
        
    if(is.null(bsl$predRec)){
      predNo <- 1
    } else{
      predNo <- max(bsl$predRec$predNo) + 1
    }
    toAdd <- data.frame(predGID, predNo, predict)
    colnames(toAdd) <- c("predGID", "predNo", "predict")
    bsl$predRec <- rbind(bsl$predRec, toAdd)

    bsl$selCriterion <- list(popID=popID, criterion="pred")
    if (exists("totalCost", bsl)) bsl$totalCost <- bsl$totalCost + bsl$costs$predCost
    return(bsl)
  }#END predict.func
  
  with(sEnv, {
    if (is.null(sharingInfo)) sharingInfo <- sims[[1]]$selCriterion$sharing
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, predictValue.func, popID=popID, trainingPopID=trainingPopID, locations=locations, years=years, sharingInfo=sharingInfo)
      sfStop()
    }else{
      sims <- lapply(sims, predictValue.func, popID=popID, trainingPopID=trainingPopID, locations=locations, years=years, sharingInfo=sharingInfo)
    }
  })
}
