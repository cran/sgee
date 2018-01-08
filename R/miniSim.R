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
#################################################################################' Mini Simulator
#' 
#' Function to do a miniature simulation testing the various grouped stagewise
#' techniques.
#' 
#' Function is currently under development and is thus not exported. Function
#' is intended to allow for simple simulation implementation to explore
#' the provided stagewise regression functions.
#' 
#' @param sampleSize Number of clusters to be used in each replicate.
#' @param clusterSize Size of each cluster.
#' @param responseCor Correlation parameter for response.
#' @param groupSparsity String specifying within group sparsity level;
#' currently takes values 'No Sparsity', 'Mod Sparsity', or
#' 'High Sparsity'.
#' @param xVariance Scalar value indicating marginal variance of the
#' covariates.
#' @param covariateCor Correlation parameter used for within group correlation.
#' @param intercept Scalar value indicating the true intercept value.
#' @param SNR Optional scalar value that allows fixing the signal
#' to noise ratio as defined as the ratio of the variance in the
#' linear predictor to the variance of the response conditioned on the
#' covariates.
#' @param techniques Vector of all techniques to be used for comparison;
#' default value of \code{c("SEE", "BiSEE", "GSEE", "HiSEE")} currently contains
#' all possible values.
#' @param stepSize Stepsize to be used in stagewise techniques.
#' @param maxNumSteps Maximum number of steps to be taken by stagewise
#' techniques.
#' @param reps Number of replicates to be used for simulation.
#' @return Simulation results.
#' @note Currently still under development, so not exported.
#' @author Gregory Vaughan
#' @keywords internal
miniSim <- function(sampleSize = 50,
                    clusterSize = 4,
                    responseCor = .4,
                    groupSparsity = "No Sparsity",
                    xVariance = 1,
                    covariateCor = .4,
                    intercept = 1,
                    SNR = 10,
                    techniques = c("SEE", "BiSEE", "GSEE", "HiSEE"),
                    stepSize = .125,
                    maxNumSteps = 400,
                    reps = 200){
    ######################
    ## Set up Coefficients
    ######################
    groupSize <- 6 ##assumed to be a multiple of 6
    numGroups <- 3
    coefs <- rep(1, groupSize)
    if(groupSparsity == "No Sparsity"){
        beta <- c(coefs,
                  rep(0, 2*groupSize))
    } else if(groupSparsity == "Mod Sparsity"){
        beta <- c(coefs[1:(groupSize/2)], rep(0, groupSize/2),
                  coefs[(groupSize/2 + 1): groupSize], rep(0, groupSize/2),
                  rep(0, groupSize))
    }else if(groupSparsity == "High Sparsity"){
        beta <- c(coefs[1:(groupSize/3)], rep(0, 2*groupSize/3),
                  coefs[(groupSize/3 + 1):(2*groupSize/3)],
                  rep(0, 2*groupSize/3),
                  coefs[(2*groupSize/3 + 1):groupSize],
                  rep(0, 2*groupSize/3))
    }

    ################
    ## Generate Data
    ################

    ## modelData is the data used to build/select the model
    modelData <- genData(numClusters = sampleSize,
                         clusterSize = clusterSize,
                         clusterRho = responseCor,
                         xVariance = xVariance,
                         numGroups = numGroups,
                         groupSize = groupSize,
                         groupRho = covariateCor,
                         beta = beta,
                         SNR = SNR,
                         intercept = intercept)
    
    
    ## predictionData is the data used to evaluate the techniques
    predictionData <- genData(numClusters = 5*sampleSize,
                              clusterSize = clusterSize,
                              clusterRho = responseCor,
                              xVariance = xVariance,
                              numGroups = numGroups,
                              groupSize = groupSize,
                              groupRho = covariateCor,
                              beta = beta,
                              SNR = SNR,
                              intercept = intercept)


    ###############
    ## Build Models
    ###############

    modelingFamily <- stats::gaussian() ## currently simulation is only for gaussian
    stoppingThreshold <- min(length(beta) ,
                             sampleSize*clusterSize-1)
    ## Make results holder
    resultArray <- array(0, c(5, length(techniques), reps))

    for(technique in techniques){

        print(technique)
        if(technique == "SEE"){
            stagewiseResult <- see(modelData$y,
                                   modelData$x,
                                   family = modelingFamily,
                                   clusterID = modelData$clusterID,
                                   maxIt = maxNumSteps,
                                   epsilon = stepSize,
                                   corstr = "exchangeable",
                                   stoppingThreshold = stoppingThreshold)
        } else if(technique == "BiSEE"){

            stagewiseResult <- bisee(modelData$y,
                                       modelData$x,
                                       family = modelingFamily,
                                       clusterID = modelData$clusterID,
                                       groupID = modelData$groupID,
                                       maxIt = maxNumSteps,
                                       epsilon = stepSize,
                                       corstr = "exchangeable",
                                       stoppingThreshold = stoppingThreshold,
                                       lambda1 = .5,
                                       lambda2 = .5)

        } else if(technique == "GSEE"){

            stagewiseResult <- gsee(modelData$y,
                                      modelData$x,
                                      family = modelingFamily,
                                      clusterID = modelData$clusterID,
                                      groupID = modelData$groupID,
                                      corstr = "exchangeable",
                                      maxIt = maxNumSteps,
                                      epsilon = stepSize,
                                      stoppingThreshold = stoppingThreshold)



    ##GSEE(modelData$Y,
    ##                                modelData$X,
    ##                                family = modelingFamily,
    ##                                clusterID = modelData$clusterID,
    ##                                groupID = modelData$groupID,
    ##                                corstr = "exchangeable",
    ##                                maxIt = maxNumSteps,
    ##                                epsilon = stepSize,
    ##                                stoppingThreshold = stoppingThreshold)

        }  else if(technique == "HiSEE"){
            stagewiseResult <- hisee(modelData$y,
                                    modelData$x,
                                    family = modelingFamily,
                                    clusterID = modelData$clusterID,
                                    groupID = modelData$groupID,
                                    corstr = "exchangeable",
                                    maxIt = maxNumSteps,
                                    epsilon = stepSize,
                                    stoppingThreshold = stoppingThreshold)
            
        }



    summary(stagewiseResult)        
    }
    

}    
    

    

