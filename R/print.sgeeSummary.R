
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
#################################################################################' \code{print} function for sgee summaries
#' 
#' Provides implementation of \code{print} function for summaries of sgee
#' objects.
#' 
#' @param x An object of the \code{sgeeSummary} class, produced by
#' applying the \code{summary} function to an object of class \code{sgee}.
#' @param ... Not currently used
#' @examples
#' 
#' #####################
#' ## Generate test data
#' #####################
#'
#' ## Initialize covariate values
#' p <- 50 
#' beta <- c(rep(2,5),
#'           c(1, 0, 1.5, 0, .5),
#'           rep(0.5,5),
#'           rep(0,p-15))
#' groupSize <- 5
#' numGroups <- length(beta)/groupSize
#' 
#' 
#' generatedData <- genData(numClusters = 50,
#'                          clusterSize = 4,
#'                          clusterRho = 0.6,
#'                          clusterCorstr = "exchangeable",
#'                          yVariance = 1,
#'                          xVariance = 1,
#'                          numGroups = numGroups,
#'                          groupSize = groupSize,
#'                          groupRho = 0.3,
#'                          beta = beta,
#'                          family = gaussian(),
#'                          intercept = 1)
#' 
#' genDF <- data.frame(generatedData$y, generatedData$x)
#' names(genDF) <- c("Y", paste0("Cov", 1:p))
#' coefMat <- hisee(formula(genDF), data = genDF,
#'                  family = gaussian(),
#'                  clusterID = generatedData$clusterID,
#'                  groupID = generatedData$groupID, 
#'                  corstr="exchangeable", 
#'                  maxIt = 50,
#'                  epsilon = .5)
#'
#' sgeeSum <- summary(coefMat)
#' print(sgeeSum)
#' 
#' @author Gregory Vaughan
#' @export 
print.sgeeSummary <- function(x, ...){
    cat("Call:\n")
    print(x$call)
    cat("\nLowest Predictive Error: ")
    cat(x$minMeasure)
    cat("\nAchieved at index: ")
    cat(x$bestIndex)    
    cat("\nWith corresponding Coefficients:\n")
    covNames <- attr(x$sgee$x, 'dimnames')[[2]]
    if(is.null(covNames)){
        covNames <- paste0("X", 1:ncol(x$sgee$x))
    }
    if(x$intercept){
        covNames <- c("(intercept)",covNames)
    }
    estimates <- matrix(x$sgee$path[x$bestIndex,],
                                 ncol = 1)
    dimnames(estimates) <- list(covNames,
                                "Estimate")
    print(estimates)   
        
    
        
}
