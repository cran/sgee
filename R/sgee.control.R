################################################################################
##
##   R package sgee by Gregory Vaughan, Kun Chen, and Jun Yan
##   Copyright (C) 2017
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
#' 
#' @return A list containing all of the parameter values.
#' 
#' @author Gregory Vaughan
#' @export
sgee.control <- function(maxIt = 200 , epsilon = 0.05,
                         stoppingThreshold =  NULL,
                         undoThreshold = 0.005,
                         interceptLimit = NULL){

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

    list(maxIt = maxIt,
         epsilon = epsilon,
         stoppingThreshold = stoppingThreshold,
         undoThreshold = undoThreshold)
}


