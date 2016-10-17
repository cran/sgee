
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
#################################################################################' Coefficient Traceplot Function
#' 
#' Function to produce the coefficent traceplot, with capabilities to
#' account for covariate groups. Used in place of the \code{plot} function.
#' 
#' plot.sgee is meant to allow for easy visualization of paths of stagewise
#' (or regularized) coefficient estimates. A great deal of flexibility is
#' provided in terms of how the plot is presented. The poenaltyFun paramter
#' allows for a penalty function to be provided (such as the $L_1$ norm)
#' to plot the coefficietn estimates against. 
#' When given the trueBeta parameter, the plot marks the paths of coefficient
#' estimates that are falsely identified as being non zero. Finally, a switch
#' for black and white versus color plots is provided (\code{color}).
#' 
#' @param x Path of coefficient Estimates.
#' @param y Optional parameter inherited from \code{plot(x,y,...)}; not used
#' with sgee.
#' @param penaltyFun Optional function that when provided results ina plot
#' of the coefficient estimates verus the corresponding penalty value.
#' When no \code{penaltyFun} value is given,
#' the plot generated is of the coefficent
#' estimates versus the iteration number.
#' @param main Optional title of plot.
#' @param xlab Label of x axis; default value is 'Iterations'.
#' @param ylab Label of y axis; default value is the beta symbol.
#' @param dropIntercept Logical parameter indicating whether the intercept
#' estimates should be dropped from the plot (i.e. not plotted). The default is
#' FALSE.
#' @param trueBeta The true coefficient values. If the true coefficient
#' values can be provided, then coefficient estimates that are false positive
#' identifications as non-zero are marked in the plot.
#' @param color Logical parameter indicating that a plot using colors
#' to differentiate coefficients is desired.
#' @param manualLinecolors Vector of desired line colors; must match dimension
#' of line colors needed (i.e. same number of colors as there are groups if
#' grouped covariates are sharing a color).
#' @param pointSpacing Space between marks used to indicate a coefficient
#' is a false positive. Spacing is measured in terms of number of indices
#' of the path matrix between marks.
#' @param ... Not currently used.
#' 
#' @note Function is intended to give a visual representation of the
#' coefficient estimates. Which x values to compare the estimates to can
#' depend on the situation, but typically the most versatile measure
#' to use is the sum of absolute values, the $L_1$ norm; especially when
#' comparing different coefficient paths from different techniques.
#' @author Gregory Vaughan
#' @examples
#' 
#' 
#' #####################
#' ## Generate test data
#' #####################
#' 
#' ## Initialize covariate values
#' p <- 50 
#' beta <- c(rep(2.4,5),
#'           c(1.3, 0, 1.7, 0, .5),
#'           rep(0.5,5),
#'           rep(0,p-15))
#' groupSize <- 1
#' numGroups <- length(beta)/groupSize
#' 
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
#'                          intercept = 0)
#' 
#' genDF <- data.frame(generatedData$y, generatedData$x)
#' coefMat <- bisee(formula(genDF),
#'                  data = genDF,
#'                  lambda1 = 0,         ##effectively see
#'                  lambda2 = 1,
#'                  family = gaussian(),
#'                  clusterID = generatedData$clusterID, 
#'                  corstr="exchangeable", 
#'                  maxIt = 200,
#'                  epsilon = .1)
#' ############################
#' ## Various options for plots
#' ############################
#' 
#' par(mfrow = c(2,2))
#'
#' ## plain useage
#' plot(coefMat, main = "Plain Usage")
#'
#' ## With penalty
#' plot(coefMat, penaltyFun = function(x){sum(abs(x))}, xlab
#' = expression(abs(abs(beta))[1]), main = "With Penalty")
#'
#' ## using true beta value to highlight misclassifications
#' plot(coefMat, trueBeta = beta, main = "ID Missclassification")
#'
#' ## black and white option
#' plot(coefMat, trueBeta = beta, color = FALSE, main =
#' "Black and White", pointSpacing = 5)
#' 
#' 
#' 
#' @export
#' @name plot.sgee
plot.sgee <- function(x,
                      y,
                      penaltyFun = NULL,
                      main = NULL,
                      xlab = "Iterations",
                      ylab = expression(beta),
                      dropIntercept = FALSE,
                      trueBeta = NULL,
                      color = TRUE,
                      manualLinecolors = NULL,
                      pointSpacing = 3,
                      ...){

    
    path <- x$path
    groupID <- x$groupID

    ## check to make sure there is even an intercept to drop
    interceptIncluded <- TRUE
    if((length(groupID) == ncol(path))){
        interceptIncluded <- FALSE
        dropIntercept <- FALSE
    }
    
    plotPoints <- rep(FALSE, ncol(path) - dropIntercept)
    pchValues <- rep(4, length(plotPoints))
    linewidths <- rep(1, ncol(path)  - dropIntercept)
    if(is.null(groupID)){
        linecolors <- grDevices::rainbow(ncol(path) - dropIntercept)
    } else{
        if(color){
            if(dropIntercept){
                linecolors <- grDevices::rainbow(length(unique(groupID)))[ groupID]
            } else {
                linecolors <- grDevices::rainbow(length(unique(groupID)) + interceptIncluded)[c(interceptIncluded, (groupID+interceptIncluded))]
            }
        } else{
            linecolors <- rep("black", ncol(path) - dropIntercept)
        }
    }
    if(color){
        if(is.null(trueBeta)){
            linetypes <- rep(1, ncol(path) - dropIntercept)
        } else{
            linetypes <- (1:2)[c(!dropIntercept, (trueBeta == 0)+1)]
        }

    } else{
        if(!is.null(trueBeta)){
                plotPoints[trueBeta ==0] <- TRUE
        }
           
        if(is.null(groupID)){
            linetypes <- rep(1, ncol(path) - dropIntercept)            
        } else{


            if(dropIntercept){
                linetypes <- (1:(length(unique(groupID))))[groupID]
            } else {
                linetypes <- (1:(length(unique(groupID))+interceptIncluded))[c(interceptIncluded, (groupID+interceptIncluded))]
            }
            
            #linetypes[linetypes  ==  3] <-4 
            #linecolors <- grDevices::rainbow(length(unique(groupID))+dropIntercept)[c(dropIntercept, (groupID+dropIntercept))]
            #cheat <- linecolors[groupID == 6]
            #linecolors[groupID == 6] <- linecolors[groupID == 2]
            #linecolors[groupID == 2] <- cheat


            ## temporarily removed, is supposed
            ## to identify important vs. unimportant groups
            ##if(!is.null(trueBeta)){
            ##    sigGroups <- unique(groupID[trueBeta !=0])
            ##    linewidths <- c(3,1)[c(dropIntercept, (!(groupID %in% sigGroups)+1))]   
            ##}

        }
           
    }

    if(!is.null(manualLinecolors)){
        if(length(manualLinecolors) == length(linecolors)){
            linecolors <- manualLinecolors
        } else{
            warning("manualLinecolors given not of appropriate length")
        }
    }


    if(is.null(penaltyFun)){
        ## Making it clear that w/o a penalty iterations
        ## are used for the x axis 
        graphics::plot(NA,NA,
                       main = main,
                       type = "n",
                       xlim = c(0,nrow(path)*1.1),
                       ylim = range(path),
                       xlab = xlab,
                       ylab = ylab)
        graphics::abline(h=0)
        ## dumb used to prevent an empty list being returned
        dumb <- sapply((1 + dropIntercept):ncol(path),function(j) { graphics::lines(1:nrow(path),
                                                                  path[,j],
                                                                  col = linecolors[j - dropIntercept],
                                                                  lty = linetypes[j - dropIntercept],
                                                                  lwd = linewidths[j - dropIntercept])
                                          
                                          if(plotPoints[j - dropIntercept]){
                                              graphics::points((1:nrow(path))[path[,j]!=0][c(TRUE, rep(FALSE, pointSpacing))], #TFF used to not overload the plot with symbols
                                                     path[,j][path[,j]!=0][c(TRUE, rep(FALSE, pointSpacing))],
                                                     col = linecolors[j - dropIntercept],
                                                     pch = pchValues[j - dropIntercept],
                                                     cex = .5)            
                                          }
                                      })

        ## If given a penaltyFun
    } else {
        if(dropIntercept){
            penaltyValues <- apply(path[,-1],
                                   MARGIN = c(1),
                                   FUN = penaltyFun)
        } else{
            penaltyValues <- apply(path,
                                   MARGIN = c(1),
                                   FUN = penaltyFun)
        }

        graphics::plot(NA,NA,
                       main = main,
                       type = "n",
                       xlim = c(0,max(penaltyValues)*1.1),
                       ylim = range(path),
                       xlab = xlab,
                       ylab = ylab,
                       ann = !(is.null(xlab) & is.null(ylab)))
        graphics::abline(h=0)
        ## dumb used to prevent an empty list being returned
        dumb <- sapply((1 + dropIntercept):ncol(path),function(j) { graphics::lines(penaltyValues,
                                                                  path[,j],
                                                                  col = linecolors[j - dropIntercept],
                                                                  lty = linetypes[j - dropIntercept],
                                                                  lwd = linewidths[j])
                                          if(plotPoints[j - dropIntercept]){
                                              graphics::points(penaltyValues[path[,j]!=0][c(TRUE, rep(FALSE, pointSpacing))], #TFFF used to not overload the plot with symbols
                                                     path[,j][path[,j]!=0][c(TRUE, rep(FALSE, pointSpacing))],
                                                     col = linecolors[j- dropIntercept],
                                                     pch = pchValues[j - dropIntercept],
                                                     cex = .5)            
                                          }
                                      })
    }
    
    
}


#' @export
#' @rdname plot.sgee
plot.sgeeSummary <- function(x,y,
                             ...){
    plot.sgee(x$sgee, y, ...)
}

