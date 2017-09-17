################################################################################
##
##   R package sgee by Gregory Vaughan, Kun Chen, and Jun Yan
##   Copyright (C) 2016-2017
##
##   This file is part of the R package sgee.
##
##   The R package sgee is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package sgee is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
#################################################################################' @export
#' @rdname bisee
gsee <- function(y, x, family,
                 clusterID, waves = NULL, groupID = 1:ncol(x),
                 corstr="independence", alpha = NULL,
                 offset = 0, 
                 intercept = TRUE,
                 control = sgee.control(maxIt = 200, epsilon = 0.05, 
                     stoppingThreshold =  min(length(y), ncol(x))-intercept,
                     undoThreshold = 0.005),
                 standardize = TRUE,
                 verbose = FALSE,
                 ...){
    
    #######################
    ## Preliminaries/set up
    #######################

    maxIt <- control$maxIt
    epsilon <- control$epsilon
    undoThreshold <- control$undoThreshold
    interceptLimit <- control$interceptLimit
    ##If the undoThreshold is >= epsilon
    ## the check wil always trigger;
    ## so this check is added to prevent
    ## an infinte loop
    if(undoThreshold >= epsilon){
        if(verbose){
            cat(paste0("****** undoThreshold too large! reducing threshold now **********\n"))
        }
        undoThreshold  <- epsilon/10
    }
    
    if(is.null(control$stoppingThreshold)){
        stoppingThreshold <- min(length(y), ncol(x))-intercept
    } else {
        stoppingThreshold <- control$stoppingThreshold
    }
        
    p <- ncol(x)

    if (is.character(family)){ 
        family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)){ 
        family <- family()
    }
    if(standardize){
        unstandardizedX <- x
        if(intercept){
            x <- scale(x)
        } else{
            x <- scale(x, center = FALSE)
        }
        
    }

    if(is.null(waves)){
        clusz <- unlist(sapply(unique(clusterID), function(x) {sum(clusterID == x)}))
        waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
    }
    
    ## currently assuming only intercept in estimation of
    ## of correlation and dispersion
    r <- 1 # number of covariates in dispersion modeling
    q <- 1 # number of covariates in correlation modeling

    ## Initial estimates for all parameters
    beta <- rep(0,p)
    phi <- stats::sd(y)^2    
    ## current initial estimation for correlation parameter
    if(is.null(alpha)){
        alpha <- 0
    }
    ## Intercept value
    beta0 <- 0
    
    meanLink <- family$linkfun
    meanLinkInv <- family$linkinv
    varianceLink <- family$variance
    mu.eta <- family$mu.eta



    ## Paths of parameter estimates
    path <- matrix(rep(0,(maxIt)*(p + intercept)), nrow = maxIt)
    ## currently assuming only intercept in estimating dispersion
    phiPath <- matrix(rep(0,(maxIt)*r), nrow = maxIt)
    ## currently assuming only intercept in estimating correlation
    alphaPath <- matrix(rep(0,(maxIt)*q), nrow = maxIt)


    numClusters <- max(clusterID)
    maxClusterSize <- max(table(clusterID))    
    mu <- rep(0,length(y))

    ##stoppedOn added to keep track of when the algorithm stops
    ## it assumes it goes the whole lenght unless stopped prematurely
    stoppedOn <- maxIt
    
    ## Working correlation matrix
    R <- genCorMat(corstr = corstr, rho = alpha, maxClusterSize = maxClusterSize)

    RInv <- solve(R)

    ##################
    ## Main Algorithim
    ##################
    cat("\n")
    oldDelta <- rep(0, length(beta)) 
    it <- 0
    while (it <maxIt){
        it <- it +1
        if(verbose){
            cat(paste0("****** Beginning iteration # ", it, " **********\n"))
        }
        GEEValues <- evaluateGEE(y = y,
                                 x = x,
                                 beta = beta,
                                 beta0 = beta0,
                                 intercept,
                                 phi = phi,
                                 offset = offset,
                                 RInv = RInv,
                                 numClusters = numClusters,
                                 clusterID = clusterID,
                                 waves = waves,
                                 meanLinkInv = meanLinkInv,
                                 mu.eta = mu.eta,
                                 varianceLink = varianceLink,
                                 corstr = corstr,
                                 maxClusterSize = maxClusterSize,
                                 interceptLimit = interceptLimit)
        ## Update Estimates
        beta0 <- GEEValues$beta0
        phi <- GEEValues$phiHat
        alpha <- GEEValues$rhoHat
        RInv <- GEEValues$RInv
        
        ## Current Values of the estimating equations
        sumMean <- GEEValues$sumMean

        ###################
        ## Update Selection
        ###################
        a <- rep(0, length(unique(groupID)))

        for (j in unique(groupID)){
          currentGroupIndex <- groupID == j
          
          aCurrent <- sumMean[currentGroupIndex]
          a[j] <- sqrt(sum(aCurrent^2))/sqrt(length(aCurrent))
        }
        ## Identify optimal update        
        delta <- which(a== max(a))
        deltaGroupIndex <- groupID == delta

        ## Check if the update is effectively undoing the last one
        if (sum(abs(oldDelta[deltaGroupIndex] + epsilon * sumMean[deltaGroupIndex]/ (a[delta]*length(sumMean[deltaGroupIndex]))))<= undoThreshold){
            if(verbose){
                cat(paste0("****** Step Undone! Reducing Stepsize **********\n"))
            }

            if(it>2){
                if (intercept){
                    beta <- path[it - 2,-1]
                } else {
                    beta <- path[it - 2,]
                }
            } else{
                beta <- rep(0, length(beta))
            }
            
            epsilon <- epsilon/2
            it <- it - 2
            oldDelta <- rep(0, length(beta))

            ##If the undoThreshold is >= epsilon
            ## the check wil always trigger;
            ## so this check is added to prevent
            ## an infinte loop
            if(undoThreshold >= epsilon){
                if(verbose){
                    cat(paste0("****** undoThreshold too large! reducing threshold now **********\n"))
                }
                undoThreshold  <- epsilon/10
            }
            
            ## If the check is passed and the update
            ## is sufficiently different from the previous
            ## update
        } else {
            oldDelta <- rep(0, length(beta))
            oldDelta[deltaGroupIndex] <- epsilon * sumMean[deltaGroupIndex]/ (a[delta]*length(sumMean[deltaGroupIndex]))

        
            ## Update estimates
            ##B_j^{[t]} = B_j^{[t-1]} - epsilon * U_j/(||U_j||_2 * sqrt(p_j))
            beta[deltaGroupIndex] <- beta[deltaGroupIndex] +  epsilon * sumMean[deltaGroupIndex]/ (a[delta]*length(sumMean[deltaGroupIndex]))

            ## Update Paths
            if(intercept){
                path[it,] <- c(beta0, beta)
            }
            else{
                path[it,] <-  beta
            }
            phiPath[it,] <- phi
            alphaPath[it,] <- alpha
            
            ###########
            ## stopping mechanism when the alogrithim has reached saturation
            if((sum(beta != 0) >= stoppingThreshold) & (it< maxIt) ){
                print("stopped on")
                print(it)
                print(a[delta])
                path[((it+1):maxIt),] <- matrix(rep(path[it,],(maxIt - it)),
                                                nrow = (maxIt - it),
                                                byrow = TRUE)
                phiPath[((it+1):maxIt),] <- matrix(rep(phiPath[it,],(maxIt - it)),
                                                   nrow = (maxIt - it),
                                                   byrow = TRUE)
                alphaPath[((it+1):maxIt),] <- matrix(rep(alphaPath[it,],(maxIt - it)),
                                                     nrow = (maxIt - it),
                                                byrow = TRUE)
                
                stoppedOn <- it
                break
            }
            
        }
    }

    

    result <- list(path = path,
                   gammaPath = phiPath,
                   alphaPath = alphaPath,
                   stoppedOn = stoppedOn,
                   maxIt = maxIt,
                   y = y,
                   x = x,
                   intercept = intercept,
                   clusterID = clusterID,
                   groupID = groupID,
                   family = family,
                   offset = offset,
                   epsilon)

    if(standardize){
        result$x <- unstandardizedX
        temp <- path
        
        if(intercept){
            temp[,-1] <- t(t(path[,-1]) /attr(x, "scaled:scale"))
            temp[,1] <- path[,1] - crossprod(t(path[,-1]),
                                             (attr(x, "scaled:center")/attr(x, "scaled:scale")))  
        } else{
            temp <- t(t(path) /attr(x, "scaled:scale"))
        }
        
        result$standardizedPath <-  path
        result$path <- temp
        result$standardizedX <- x
    }

    class(result) <-  "sgee"

    result
}
