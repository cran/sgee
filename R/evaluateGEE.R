################################################################################
##
##   R package sgee by Gregory Vaughan, Kun Chen, and Jun Yan
##   Copyright (C) 2016-2018
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
################################################################################## Function intended to support gBoost, sgBoost, and hBoost
## functions by evaluating the GEEs 
evaluateGEE <- function(y,
                        x,
                        beta,
                        beta0,
                        intercept,
                        phi,
                        offset = 0,
                        RInv,
                        numClusters,
                        clusterID,
                        waves = NULL,
                        meanLink = NULL,
                        meanLinkInv,
                        mu.eta,
                        varianceLink,
                        corstr,
                        maxClusterSize,
                        scaleObjective = FALSE,
                        interceptLimit){

    p <- sum(beta!=0) +  intercept
    sumMean <- rep(0,length(beta))
    ##HMean <- matrix(0, ncol = length(beta), nrow = length(beta))
    



    
    
    ##################
    ## update intercept
    ##################
    if(intercept){
        interceptEstimatingEquation <- function(beta0){
            eta <- x %*% beta + beta0 +offset
            mu <- meanLinkInv(eta)
            Variancei <- varianceLink(mu) * phi

            sumIntercept <- 0
            ##HIntercept <- 0
            
            for(i in unique(clusterID)){
                currentClusterIndex <- clusterID == i
                
                ## pearsonResidual = viHalfInv %*% (Y-mu)
                pearsonResiduals <- (y[currentClusterIndex] - mu[currentClusterIndex]) /sqrt(Variancei[currentClusterIndex])

            
                ## VhalfDiMean = V^(-1/2) %*% DiMean 
                VhalfDiIntercept <- matrix(mu.eta(eta[currentClusterIndex]) * (1/sqrt(Variancei[currentClusterIndex])), ncol = 1)
                
                ## DPlus = R^(-1) %*% VhalfDiMean
                ##if(is.null(waves)){
                ##     DPlus <- crossprod(RInv,VhalfDiIntercept)
                ## } else{
                    curWave <- waves[clusterID == i]
                    DPlus <- crossprod(RInv[curWave, curWave],VhalfDiIntercept)
                ##}
                
            
                ## sum is used since DPlus and VhalfDiIntercept will be vectors here
                
                sumIntercept <- sumIntercept + sum(DPlus* pearsonResiduals)
                ##HIntercept <- HIntercept + sum(VhalfDiIntercept * DPlus)
                
            }
            
            return(sumIntercept)
        }

        ## update intercept Term
        if(is.null(interceptLimit)){
            ## Uniroot allows for extending the interval being investigated
            ## but this extending functionality is new, and so there may be
            ## errors in older versions of R.
            ## set interceptLimit to a positive number in the sgeeControl to
            ## avoid this issue

            beta0 <- stats::uniroot(interceptEstimatingEquation,
                                    c(-10, 10),
                                    ##extendInt causes error in older R versions
                                    #trace = 2,
                                    extendInt = "yes" 
                                    )$root            
        } else {
            beta0 <- stats::optimize(function(y){(interceptEstimatingEquation(y))^2},
                                     c(-interceptLimit, interceptLimit))$minimum
        }
        
    }

    #########################
    ## Estimate phi and alpha
    #########################
    eta <- x %*% beta + beta0 + offset
    mu <- meanLinkInv(eta)
    Variancei <- varianceLink(mu) * phi 

    phiNumerator <- phiDenominator <- 0
    if(corstr == "exchangeable"){
        ##numObservations <-0
        rhoNumerator <- rhoDenominator <- 0
    } else if(corstr == "ar1"){
        rhoNumerator <- rhoDenominator <- 0
        ## experimental approach the was not completed
        ## numPairs <- 0
        ## for (i in unique(clusterID)){
        ##     currentClusterIndex <- clusterID == i 
        ##     ni <- sum(currentClusterIndex)
        ##     if(ni >1){
        ##         ## ncol(combn(ni,2)) = n choose 2
        ##         numPairs <- numPairs + ncol(combn(ni,2))
        ##     }
        ## }
        
        ## pearsonPairs <- rep(0, numPairs)
        ## ##|t-t'| values, where t and t' are in-cluster indecies
        ## lagVals <- rep(0, numPairs)

        ## ## used to keep track of position in pearsonPairs and lagVals
        ## pairsIndex <- 1
        
        
    }
    
    for (i in unique(clusterID)){
        currentClusterIndex <- clusterID == i 
        
        
        ## pearsonResidual = viHalfInv %*% (Y-mu)
        pearsonResiduals <- (y[currentClusterIndex] - mu[currentClusterIndex]) /sqrt(Variancei[currentClusterIndex]/phi)
        
        ##numerator for Phi Estimator
        phiNumerator <- phiNumerator + sum(pearsonResiduals^2)

        ## product of different pearson residual pairs
        z <- tcrossprod(pearsonResiduals, pearsonResiduals)

        ## number of observation in this cluster
        ni <- sum(currentClusterIndex)
        if(corstr == "exchangeable"){
            ##numerator for correlation (rho) estimator
            rhoNumerator <- rhoNumerator + sum(z[upper.tri(z)])
            
            ## denominator for correlation estimator
            ##numObservations <- numObservations + ni
            rhoDenominator <- rhoDenominator + 0.5*(ni-1)*(ni)
        } else if(corstr =="ar1"){
            ##numerator for correlation (rho) estimator
            rhoNumerator <- rhoNumerator + sum(diag(z[(1:(ni-1)),(2:ni)]))
            
            ## denominator for correlation estimator
            ##numObservations <- numObservations + ni
            rhoDenominator <- rhoDenominator + (ni-1)
         
            ## experimental approach the was not completed
            ## curWave <- waves[clusterID == i]
            
            ## pearsonPairs[pairsIndex:(pairsIndex + ncol(combn(ni,2)) -1)] <- z[upper.tri(z)]

            ## lagMatrix <- abs(outer(curWave, curWave, "-"))
            ## ##|t-t'| values, where t and t' are in-cluster indecies
            ## lagVals[pairsIndex:(pairsIndex + ncol(combn(ni,2)) -1)] <- lagMatrix[upper.tri(lagMatrix)]


            ## ## ncol(combn(ni,2)) = n choose 2
            ## pairsIndex <- pairsIndex + ncol(combn(ni,2))
        }
        
    }

    phiDenominator <- (length(clusterID) - p)
    phi <- phiNumerator/phiDenominator

    if(corstr == "independence"){
        rho <- 0
    }else if(corstr == "exchangeable"){
        rhoDenominator <- (rhoDenominator - p) *(phiNumerator/phiDenominator)
        rho <- rhoNumerator/rhoDenominator
    }else if(corstr == "ar1"){
        rhoDenominator <- (rhoDenominator - p) *(phiNumerator/phiDenominator)
        rho <- rhoNumerator/rhoDenominator
    }
    
    ## a catch to ensure that rho is <1, which can occassionally happen
    ## in the beginning due to poor initial estimates
    if(rho >= 1){
        warning("Rho had to be forced to be less than 1")
        rho <- 0.99
    }
    R <- genCorMat(corstr = corstr, rho = rho, maxClusterSize = maxClusterSize)
    RInv <- solve(R)
    
    ###########################
    ## Estimate mean parameters
    ###########################
    eta <- x %*% beta + beta0 + offset
    mu <- meanLinkInv(eta)
    Variancei <- varianceLink(mu) * phi 


    for (i in unique(clusterID)){       
        currentClusterIndex <- clusterID == i 
        ## mean paramter estimating equation
        
        ## pearsonResidual = viHalfInv %*% (Y-mu)
        pearsonResiduals <- (y[currentClusterIndex] - mu[currentClusterIndex]) /sqrt(Variancei[currentClusterIndex]/phi)
        
        
        ## VhalfDiMean = V^(-1/2) %*% DiMean 
        VhalfDiMean <- mu.eta(eta[currentClusterIndex]) * x[currentClusterIndex,, drop = FALSE] * (1/sqrt(Variancei[currentClusterIndex]))
        
        ## DPlus = R^(-1) %*% VhalfDiMean
        ## if(is.null(waves)){
        ##    DPlus <- crossprod(RInv,VhalfDiMean)
        ## } else{
            curWave <- waves[clusterID == i]
            DPlus <- crossprod(RInv[curWave, curWave],VhalfDiMean)
        ## }
        
        ## SumMean = sum(t(DiMean)V^(-1/2)R^(-1)V^(-1/2)(Y-mu))
        ## HMean = sum(t(DiMean)V^(-1/2)R^(-1)V^(-1/2)DiMean)
        sumMean <- sumMean + crossprod(DPlus, pearsonResiduals)
        ##HMean <- HMean + crossprod(VhalfDiMean,DPlus)        
    }
      
    return(list(beta0 = beta0,
                sumMean = sumMean,
                ##HMean = HMean,
                phiHat = phi, ##phiNumerator/phiDenominator
                rhoHat = rho, ##rhoNumerator/rhoDenominator
                RInv = RInv))
    
}
