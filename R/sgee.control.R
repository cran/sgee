################################################################################
##
##   R package sgee by Gregory Vaughan, Kun Chen, and Jun Yan
##   Copyright (C) 2017-2018
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
#' Auxiliary for Controlling SGEE fitting
#' 
#' Auxiliary function for \code{sgee} fitting functions. Specifies
#' parameters used by all \code{sgee} fitting functions in terms
#' of the path generation; i.e. step size \code{epsilon}, maximum number of
#' iterations \code{maxIt}, and the threshold for premature stopping
#' \code{stoppingthreshold}.
#' 
#' 
#' @param epsilon Step size to be used when incrementing coefficient value(s)
#' in each iteration. Default is 0.05.
#' @param maxIt Maximum number of iterations of the stagewise algorithm to be
#' executed. Default is 200.
#' @param stoppingThreshold An integer value that indicates the maximum number
#' of allowed covariates in the model. Once the algorithm has reached the value
#' of \code{stoppingThreshold}, the algorithm will stop without completing any
#' remaining iterations. The number of covariates to be included cannot exceed
#' the number of observations. The default value is typically
#' the minimum of the number
#' of covariates and the number of observations, minus 1 if an intercept is
#' included.
#' @param undoThreshold A small value used to determine if consecutive
#' steps are sufficiently different. If consecutive steps effectively undo each
#' other (as indicated by having a sum with an absolute value less than
#' \code{undoThreshold}), then the steps are repeated and the stepsize is
#' reduced. A negative value for \code{undoThreshold} effectively prevents
#' this step. \code{undoThreshold} should only be big enough to allow for
#' some rounding error in steps and should be much smaller than the step size.
#' Default value is 0.005.
#' @param interceptLimit sgee functions make use of the extendInt parameter
#' of uniroot to estimate the intercept in each iteration. This parameter
#' was recently implemented and thus may cause issues with older versions of R.
#' If a value is given for \code{interceptLimit}, then this extendInt parameter
#' is bypassed and a solution for the intercept estimating equation is sought
#' out between negative \code{interceptLimit} and positive
#' \code{interceptLimit}. The default value of \code{NULL} uses the
#' extendInt functionality.
#' @param stochastic A numeric value between 0 (exclusive) and 1 (inclusive)
#' to indicate what proportion of the data should be subsampled
#' in the stochastic implementation of stagewise approached. The default value
#' of 1 implements the standard deterministic approach where no subsampling
#' is done.
#' @param sampleProb A user provided value dictating the
#' probability distribution for stochastic stagewise approaches.
#' \code{sampleProb} can be provided as 1) a vector of
#' fixed values of length equal to the resposne vector y, 2) a function
#' that takes in a list of values (full list of values given in details)
#' and returns a vector of length equal to the response vector y, or 3) the
#' default value of \code{NULL}, which results in a uniform distribution
#' @param reSample Parameter indicating how frequently a subsample is
#' collected in stochastic stagewise approaches. If reSample == 1 then
#' a subsample is collected every iteration,if reSample == 2 a subsample
#' is collected every two (i.e every other) iteration. If reSample == 0,
#' (the default value) then a subsample is only collected once
#' before any iterations have been done.
#' @param withReplacement a Logical value indicating if the subsampling in
#' stochastic stagewise approaches should be done with or without replacement.
#' Default values is \code{FALSE}.
#' 
#' @return A list containing all of the parameter values.
#' 
#' @author Gregory Vaughan
#' @export
sgee.control <- function(maxIt = 200 , epsilon = 0.05,
                         stoppingThreshold =  NULL,
                         undoThreshold = 0.005,
                         interceptLimit = NULL,
                         stochastic = 1,
                         sampleProb = NULL,
                         reSample = 1,
                         withReplacement = FALSE){

    ## some basic checks
    if(epsilon <= 0){
        stop("epsilon must be >0")
    }
    if(epsilon >=1){
        warning("It is advised that epsilon be sufficiently small")
    }
    if(maxIt <= 0){
        stop("maxIt must be >0")
    }
    if(!is.null(stoppingThreshold)){
        if(stoppingThreshold<1){
            stop("stoppingThreshold must be >1")
        }
    }

    if(!is.null(interceptLimit)){
        if(interceptLimit<=0){
            stop("interceptLimit must be >0")
        }
    }

    ##############
    ## stochastic value indicates whether subsampling is used,
    ## must be between 0 (exclusive) and 1 (inclusive)
    if(stochastic>1 | stochastic<=0){
        stop("paramter 'stochastic' indicates percentage of data being subsample and must be greater than 0 and less than or equal to 1")
    }
    if(reSample <0 | reSample %% 1 != 0){
        stop("parameter 'reSample' indicates how frequently a subsample is collected (i.e. reSample == 1 implies every iteration, reSample == 2 implies every other iteration). Therefore, reSample must be a non-negative integer")
    }
    

    list(maxIt = maxIt,
         epsilon = epsilon,
         stoppingThreshold = stoppingThreshold,
         undoThreshold = undoThreshold,
         stochastic = stochastic,
         sampleProb = sampleProb,
         reSample = reSample,
         withReplacement = withReplacement)
}


