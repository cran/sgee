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
#################################################################################' Correlation Matrix Generator.
#' 
#' Function that generates a correlation matrix of a predefined type  and size
#' given appropriate correlation parameter(s), \code{rho}.
#' 
#' 
#' @param corstr Structure of correlaiton matrix to be generated;
#' 'independence', 'exchangeable', 'ar1', and 'unstructured' currently
#' implemented.
#' @param rho Correlation parameter; assumed to be of length 1 or
#' \code{maxClusterSize * (maxClusterSize - 1) /2}. 
#' @param maxClusterSize size of the correlation matrix being generated.
#' @return A correlation matrix of form matching corstr and of size
#' maxClusterSize.
#' 
#' @note Mostly intended for internal use, but could be useful to user.
#' Therefore, the function is exported.
#' @author Gregory Vaughan
#' @examples
#' 
#' 
#' ## Generates Correlation Matricies easily
#' ## When corstr = "independence", the value of rho
#' ## is irrelevant
#' mat1 <- genCorMat(corstr = "independence", rho = .1, maxClusterSize = 3) 
#' 
#' ## Exchangeable
#' mat2 <- genCorMat(corstr = "exchangeable", rho = .3, maxClusterSize = 2) 
#' 
#' ## AR-1
#' mat3 <- genCorMat(corstr = "ar1", rho = .4, maxClusterSize = 4) 
#' 
#' ## unstructured
#' mat3 <- genCorMat(corstr = "unstructured",
#'                   rho = c(.3,.2,.1),
#'                   maxClusterSize = 3) 
#'
#' 
#' @export genCorMat
genCorMat <-
function(corstr="independence", rho, maxClusterSize = 0){


    corstr <- tolower(corstr)
    R <- diag(.5, maxClusterSize)
    upperTriangle <- upper.tri(R)
    
    if(maxClusterSize == 0){
        stop("maxClusterSize must be greater than 0")
    }


        
    if(corstr == "independence"){
        R <- diag(1, maxClusterSize)
    }
    else if(corstr == "exchangeable" | corstr == "ex"){
        if (length(rho) != 1){
            print(paste0("length of rho: ", length(rho)))
            stop("length of rho must 1 for corstr = exchangeable.") 
        }
        R[upperTriangle] <- rho
        R <- R + t(R)
    }
    else if(corstr == "ar1"){
        if (length(rho) != 1){
            print(paste0("length of rho: ", length(rho)))
            stop("length of rho must 1 for corstr = ar1.") 
        }
        powers <- abs(outer(1:maxClusterSize, 1:maxClusterSize, "-"))
        R[upperTriangle] <- rho
        R <- R + t(R)

        ## Absolute value of corMat entries
        absR <- abs(R)^powers

        ## code that allows for negative rho values
        signPowers <- !diag(1, maxClusterSize)
        R <- (sign(R)^signPowers) * absR
    }
    else if(corstr == "unstructured"){
        if (length(rho) != maxClusterSize * (maxClusterSize-1) /2){
            print(paste0("length of rho: ", length(rho)))
            if (length(rho) ==1){
                print("univariate rho value is equivalent to corstr = exchangeable")
            }
            stop("length of rho must either be maxClusterSize * (maxClusterSize-1) /2 for corstr = unstructured.") 
        }

        R[upperTriangle] <- rho
        R <- R + t(R)
        
    }
    else{
        stop("invalid corstr value.")
    }
    

    
    if(!all(abs(R[upperTriangle])<=1)){
        print(R)
        stop("invalid working correlation structure generated with improper values.")
    }


    return(R)
    
}
