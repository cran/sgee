
################################################################################
##
##   R package sgee by Gregory Vaughan, Kun Chen, and Jun Yan
##   Copyright (C) 2016
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
see <- function(y, x,  family,
                clusterID,
                waves = NULL,
                corstr="independence", alpha = NULL,
                intercept = TRUE,
                 offset = 0,
                control = sgee.control(maxIt = 200, epsilon = 0.05, 
                    stoppingThreshold =  min(nrow(y), ncol(x))-intercept),
                standardize = TRUE,
                ...){
    #######################
    ## Preliminaries/set up
    #######################
    
    maxIt <- control$maxIt
    epsilon <- control$epsilon
    if(is.null(control$stoppingThreshold)){
        stoppingThreshold <- min(nrow(y), ncol(x))-intercept
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
    
    ## currently assuming only intercept in estimation
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
    ## it assumes it goes the whole length unless stopped prematurely
    stoppedOn <- maxIt

    ## Working correlation matrix
    R <- genCorMat(corstr = corstr, rho = alpha, maxClusterSize = maxClusterSize)

    RInv <- solve(R)

    ##################
    ## Main Algorithim
    ##################
    for (it in 1:maxIt){
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
                                 maxClusterSize = maxClusterSize)
        ## Update Estimates
        beta0 <- GEEValues$beta0
        phi <- GEEValues$phiHat
        alpha <- GEEValues$rhoHat
        RInv <- GEEValues$RInv

        ## Current values of the estimating equations
        sumMean <- GEEValues$sumMean

        ## Identify optimal update
        a <- abs(sumMean)   
        delta <- which(a== max(a))

        ## Update estimates
        beta[delta] <- beta[delta] +epsilon*sign(sumMean[delta]) 

        ## Update the paths
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

    

    result <- list(path = path,
                   gammaPath = phiPath,
                   alphaPath = alphaPath,
                   stoppedOn = stoppedOn,
                   maxIt = maxIt,
                   y = y,
                   x = x,
                   intercept = intercept,
                   clusterID = clusterID,
                   groupID = 1:ncol(x),
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

    class(result) <- "sgee"

    result
}