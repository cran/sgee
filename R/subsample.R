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
#' subsample
#' 
#' Internal function to execute the subsampling component of
#' the stochastic stagewise approach. If a user provides a \code{stochastic}
#' value between 0 and 1, it is assumed that some proportion of subsampling
#' is desired. The \code{samplingDistCalculation} function calculates the
#' distribution of the clusters and the \code{subsample} function uses that
#' distribution to draw the actual subsample. 
#' 
#' @param sampleDist A vector whose length is equal to the number of clusters
#' that indicates the probability of sampling each cluster
#' @param sampleSize A scalar value indicating how larger of a subsample is
#' being drawn
#' @param withReplacement A logical value indicating whether the
#' subsampling is beign done with or without replacement
#' @param clusterIDs A vector of all of the UNIQUE cluster IDs
#' @param clusterID A vector of length equal to the number of observations
#' indicating which cluster each observation is in
#'
#' 
#' @return A list with two variables: \code{subSampleIndicator}, which
#' indicates which observations are in the current subsample, and
#' \code{clusterIDCurr}, which indicates the clusterID for the subsample.
#'
#' @note Internal function.
#'
#' While most of the subsample can be determined from the
#' \code{subSampleIndicator}, the \code{clusterIDCurr} value has to be
#' constructed inside the \code{subsample function} as the way the cluster
#' IDs is handled is different depending o n whether we are sampling with
#' or without replacement. 
#' 
#' 
#' @author Gregory Vaughan
#' @keywords internal
subsample <- function(sampleDist,
                      sampleSize,
                      withReplacement,
                      clusterIDs,
                      clusterID){

    IDsSample <- sample( clusterIDs, sampleSize, prob = sampleDist, replace = withReplacement)
    ## special counters need to be set up
    ## to sample clusters with replacement
    ## properly
    if(withReplacement){
        subSampleIndicator <- numeric(0)
        clusterIDCurr <- numeric(0)
        counter <- 0
        for(id in IDsSample){
            counter <- counter + 1
            subSampleIndicator <- c(subSampleIndicator, which(clusterID == id))
            clusterIDCurr <- c(clusterIDCurr, rep(counter, sum(clusterID == id))) 
        }
    } else{
        subSampleIndicator <- clusterID %in% IDsSample
        clusterIDCurr <- clusterID[subSampleIndicator]
    }


    
    subsampleOutput <- list(subSampleIndicator = subSampleIndicator,
                            clusterIDCurr = clusterIDCurr)
    subsampleOutput
}
