
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
                        estimateRho,
                        scaleObjective = FALSE){

    p <- sum(beta!=0) +  intercept
    sumMean <- rep(0,length(beta))
    HMean <- matrix(rep(0,length(beta)^2),nrow = length(beta))
    

    numObservations <-0

    
    
    ##################
    ## update intercept
    ##################
    if(intercept){
        interceptEstimatingEquation <- function(beta0){
            eta <- x %*% beta + beta0 +offset
            mu <- meanLinkInv(eta)
            Variancei <- varianceLink(mu) * phi

            sumIntercept <- 0
            HIntercept <- 0

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
                HIntercept <- HIntercept + sum(VhalfDiIntercept * DPlus)
                
            }
            
            return(sumIntercept)
        }

        ## update intercept Term
        beta0 <- stats::uniroot(interceptEstimatingEquation,
                                c(-100, 100),
                                ##extendInt causes error in older R versions 
                                extendInt = "yes" 
                                )$root
    }

    #########################
    ## Estimate phi and alpha
    #########################
    eta <- x %*% beta + beta0 + offset
    mu <- meanLinkInv(eta)
    Variancei <- varianceLink(mu) * phi 

    phiNumerator <- phiDenominator <- rhoNumerator <- rhoDenominator <-0    
    for (i in unique(clusterID)){
        currentClusterIndex <- clusterID == i 
        
        ## pearsonResidual = viHalfInv %*% (Y-mu)
        pearsonResiduals <- (y[currentClusterIndex] - mu[currentClusterIndex]) /sqrt(Variancei[currentClusterIndex]/phi)
        
        ##numerator for Phi Estimator
        phiNumerator <- phiNumerator + sum(pearsonResiduals^2)
        
        ##numerator for correlation (rho) estimator
        z <- tcrossprod(pearsonResiduals, pearsonResiduals)
        
        rhoNumerator <- rhoNumerator + sum(z[upper.tri(z)])
        
        ## denominator for correlation estimator
        ni <- sum(currentClusterIndex)
        numObservations <- numObservations + ni
        rhoDenominator <- rhoDenominator + 0.5*(ni-1)*(ni)
        
    }

    phiDenominator <- (length(clusterID) - p)
    rhoDenominator <- (rhoDenominator - p) *(phiNumerator/phiDenominator)


    phi = phiNumerator/phiDenominator
    rho = rhoNumerator/rhoDenominator    
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
        HMean <- HMean + crossprod(VhalfDiMean,DPlus)        
    }
      
    return(list(beta0 = beta0,
                sumMean = sumMean,
                HMean = HMean,
                phiHat = phiNumerator/phiDenominator,
                rhoHat = rhoNumerator/( rhoDenominator),
                RInv = RInv))
    
}
