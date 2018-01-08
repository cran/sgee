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
################################################################################
#' samplingDistCalculation
#' 
#' Internal function to set up subsampling distribution
#' to execute the stochastic version of a stagewise approach.
#' The subsampling is coducted at the cluster level, not the
#' individual observation level. Sampling probabilities are first
#' calculated or provided for each observation individually, and
#' then the sampling probability for each cluster is taken to be
#' the average probability across all observations in the cluster.
#' 
#' 
#' @param sampleProb A user provided value for the probability associated
#' with each observation. \code{sampleProb} can be provided as 1) a vector of
#' fixed values of length equal to the resposne vector y, 2) a function
#' that takes in a list of values (full list of values given in details)
#' and returns a vector of length equal to the response vector y, or 3) the
#' default value of \code{NULL}, which results in a uniform distribution
#' @param y The vector of the response values provided to the original
#' stagewise function
#' @param x The covariate matrix provided to the original stagewise function
#' @param clusterID The vector of cluster ID numbers provided to the original
#' stagewise function
#' @param waves The waves parameter identifying the order of observations
#' within the clusters that is provided to the original stagewise function
#' @param beta The vector of the current estimates of the coefficients
#' @param beta0 The current estimate of the intercept
#' @param phi Current estimate of the scale parameter
#' @param alpha Current estimate of the parameter affecting the within
#' cluster correlation
#' @param offset offset in the linear predictor provided to the original
#' stagewise function
#' @param meanLinkInv The link inverse function from the \code{family}
#' object provided to the original stagewise function indicating what family
#' of mean and variance structure is assumed
#' @param varianceLink The variance link function from the \code{family}
#' object provided to the original stagewise function indicating what family
#' of mean and variance structure is assumed
#' @param corstr The structure of the working correlation matrix that was
#' provided to the original stagewise function
#' @param mu.eta Derivative function of mu, the conditional mean of the
#' response, with respect to eta, the linear predictor, from the \code{family}
#' object provided to the original stagewise function indicating what family
#' of mean and variance structure is assumed
#'
#' 
#' @return The sampling distribution probabilities to be used for the sub
#' sampling. distribution is provided as a vector with length equal to the
#' number of clusters.
#'
#' @note Internal function.
#'
#' The function provided to \code{sampleProb} (through the
#' \code{sgee.control} function) needs to calculate
#' probabilities for each observation in the response vector \code{y}.
#' How these calculations are done is up to the user and the following
#' values are provided to the \code{sampleProb} function as a list called
#' \code{values}: \code{y}, \code{x}, \code{clusterID}, \code{waves},
#' \code{beta}, \code{beta0}, \code{phi}, \code{alpha}, \code{offset},
#' \code{meanLinkInv}, \code{varianceLink}, \code{corstr}, \code{mu.eta}.
#' additionally, all of the values produced by \code{sampleProb} need to be
#' non-negative.
#' 
#' 
#' @author Gregory Vaughan
#' @keywords internal
samplingDistCalculation <- function(sampleProb,
                                    y,
                                    x,
                                    clusterID,
                                    waves,
                                    beta,
                                    beta0,
                                    phi,
                                    alpha,
                                    offset,
                                    meanLinkInv,
                                    varianceLink,
                                    corstr,
                                    mu.eta){

    ## default subsampling distribnution is uniform
    if(is.null(sampleProb)){
      ## turns out we actually want to keep it null
      ## sample()still does a uniform sampling,
      ## but much faster
        #sampleDist <- rep(1/length(y), length(y))
        sampleDist <- NULL
    } else if(is.function(sampleProb)){
        sampleDist <- sampleProb(values = list(y = y, x = x, clusterID = clusterID, waves = waves, beta = beta, beta0 = beta0, phi = phi, alpha = alpha, offset = offset, meanLinkInv = meanLinkInv, varianceLink = varianceLink, corstr = corstr, mu.eta = mu.eta ))
        sampleDist <- sampleDist/sum(sampleDist)
        if(any(sampleDist < 0)){
            stop("sampleProb Function must be a non-negative function")
        }
    } else if(is.vector(sampleProb)){
        if(is.numeric(sampleProb) & length(sampleProb) == length(y)){ 
            sampleDist <-sampleProb
        } else {
            stop("sampleProb is a vector, but is either not numeric or has a length different than y")
        }
    } else{
        stop("sampleProb must either be a non-negative function or a numeric vector of length equal to y")
    }

    
    ## taking the average prob across the cluster to get the cluster prob
    ## if sampleDist is NULL, then the default behavior of sample()
    ## uses uniform distribution
    if(!is.null(sampleDist)){
        sampleDist <- by(sampleDist, INDICES = clusterID, FUN = mean) 
    }
    sampleDist
    
}
