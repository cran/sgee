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
#################################################################################' Bi-level Stagewise Estimating Equations Implementation
#'
#' 
#' Function to perform BiSEE, a Bi-Level Boosting / Functional
#' Gradient Descent / Forward Stagewise regression in
#' the grouped covariates setting using
#' Generalized Estimating Equations
#' 
#' Function to implement BiSEE, a stagewise regression approach
#' that is designed to perform bi-level selection in the context of
#' Generalized Estimating Equations. Given a response \code{y} and
#' a design matrix \code{x}
#' (excluding intercept) BiSEE generates a path of stagewise regression
#' estimates for each covariate based on the provided step size epsilon,
#' and tuning parameters \code{lambda1} and \code{lambda2}. When
#' \code{lambda1 == 0} or \code{lambda2 == 0}, the simplified versions
#' of \code{bisee} called \code{see} and \code{gsee}, respectively,
#' will be called.
#'
#' The resulting path can then be analyzed to determine an optimal
#' model along the path of coefficient estimates. The
#' \code{summary.sgee} function provides such
#' functionality based on various
#' possible metrics, primarily focused on the Mean Squared Error.
#' Furthermore, the \code{plot.sgee} function can be used to examine the
#' path of coefficient estimates versus the iteration number, or some
#' desired penalty.
#'
#' \code{bisee} makes use of the function uniroot in the stats package.
#' The \code{extendInt} parameter for \code{uniroot} is used, which may
#' cause issues for older versions of R.
#'
#' @param y Vector of response measures that corresponds with modeling family
#' given in 'family' parameter. \code{y} is assumed to be the same length as
#' \code{clusterID} and is assumed to be organized into clusters as dictated by
#' \code{clusterID}.
#' @param x Design matrix of dimension \code{length(y)} x nvar,
#' the number of variables, where each row is
#' represents an obersvation of predictor variables. Assumed to be scaled.
#' @param family Modeling family that describes the marginal distribution of
#' the response. Assumed to be an object such as \code{gaussian()} or
#' \code{poisson()}.
#' @param clusterID Vector of integers that identifies the clusters of response
#' measures in \code{y}. Data and \code{clusterID} are assumed to
#' 1) be of equal lengths,
#' 2) sorted so that observations of a cluster are in contiguous rows, and 3)
#' organized so that \code{clusterID} is a vector of consecutive integers.
#' @param groupID Vector of integeres that identifies the groups of the
#' covariates/coefficients (i.e. the columns of \code{x}).  \code{x}
#' and \code{groupID} are
#' assumed 1) to be of corresponding dimension, (i.e.
#' \code{ncol(x) == length(groupID)}),
#' 2) sorted so that groups of covariates are in contiguous
#' columns, and 3) organized so that \code{groupID} is a vector of consecutive
#' integers.
#' @param corstr A character string indicating the desired working correlation
#' structure. The following are implemented : "independence" (default value),
#' "exchangeable", and "ar1".
#' @param alpha An initial guess for the correlation parameter value
#' between -1 and 1 . If left NULL (the default), the initial estimate is 0.
#' @param lambda1 Mixing parameter used to indicate weight of $L_2$ Norm
#' (group selection). While not necessary, \code{lambda1} and
#' \code{lambda2} are
#' best set to sum to 1 as only the weights relative to each other matter.
#' Default value is set to .5.
#' @param lambda2 Mixing parameter used to indicate weight of $L_1$ Norm
#' (individual selection). While not necessary, \code{lambda1} and
#' \code{lambda2} are best set to sum
#' to 1 as only the weights relative to each other matter.  Default value is
#' set to \code{1-lambda1}.
#' @param intercept Binary value indicating where an intercept term is
#' to be included in the model for estimation. Default is to include an
#' intercept.
#' @param offset Vector of offset value(s) for the linear predictor.
#' \code{offset}
#' is assumed to be either of length one, or of the same length as \code{y}.
#' Default is to have no offset.
#' @param control A list of parameters used to contorl the path generation
#' process; see \code{sgee.control}.
#' @param standardize A logical parameter that indicates whether or not
#' the covariates need to be standardized before fitting.
#' If standardized before fitting, the unstandardized
#' path is returned as the default, with a \code{standardizedPath} and
#' \code{standardizedX} included
#' separately. Default value is \code{TRUE}.
#' @param verbose Logical parameter indicating whether output should be produced
#' while bisee is running. Default value is FALSE.
#' @param ... Not currently used
#' 
#' @return Object of class \code{sgee} containing the path
#' of coefficient estimates,
#' the path of scale estimates, the path of correlation parameter
#' estimates, the iteration at which BiSEE terminated, and initial regression
#' values including \code{x}, \code{y}, code{family}, \code{clusterID},
#' \code{groupID}, \code{offset}, \code{epsilon}, and \code{numIt}.
#' 
#' @note Function to execute BiSEE technique. Note that \code{lambda1}
#' and \code{lambda2} are tuning parameters. Though it is advised to
#' fix \code{lambda1 + lambda2 = 1}, this is not necessary. These parameters
#' can be tuned using various approaches including cross validation.
#' @author Gregory Vaughan
#' @references G. Vaughan, R. Aseltine, K. Chen & J. Yan (2016). Stagewise
#' Generalized Estimating Equations with Grouped Variables. Department of
#' Statistics, University of Connecticut. Technical Report 16-09.
#'
#' Wolfson, J. (2011). EEBoost: A general method for prediction
#' and variable selection based on estimating equations. Journal of the
#' American Statistical Association 106, 296--305.
#'
#' Tibshirani, R. J. (2015). A general framework for fast stagewise
#' algorithms. Journal of Machine Learning Research 16, 2543--2588.
#'
#' Simon, N., Friedman, J., Hastie, T., and Tibshirani, R. (2013). A
#' sparse-group lasso. Journal of Computational and Graphical
#' Statistics 22, 231--245.
#' @examples
#' 
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
#' 
#' ## Perform Fitting by providing y and x values
#' coefMat1 <- bisee(y = generatedData$y, x = generatedData$x,
#'                  family = gaussian(),
#'                  clusterID = generatedData$clusterID,
#'                  groupID = generatedData$groupID, 
#'                  corstr = "exchangeable", 
#'                  control = sgee.control(maxIt = 50, epsilon = 0.5),
#'                  lambda1 = .5,
#'                  lambda2 = .5,
#'                  verbose = TRUE)
#'
#'
#' ## Perform Fitting by providing formula and data
#' genDF <- data.frame(generatedData$y, generatedData$x)
#' names(genDF) <- c("Y", paste0("Cov", 1:p))
#' coefMat2 <- bisee(formula(genDF), data = genDF,
#'                  family = gaussian(),
#'                  subset = Y <1.5,
#'                  waves = rep(1:4, 50), 
#'                  clusterID = generatedData$clusterID,
#'                  groupID = generatedData$groupID, 
#'                  corstr = "exchangeable",
#'                  control = sgee.control(maxIt = 50, epsilon = 0.5),
#'                  lambda1 = 0.5,
#'                  lambda2 = 0.5,
#'                  verbose = TRUE)
#' 
#' par(mfrow = c(2,1))
#' plot(coefMat1)
#' plot(coefMat2)
#' 
#' @export bisee
#' @name bisee
NULL

#' @export
#' @rdname bisee 
bisee <- function(y, ...) UseMethod("bisee")

#' @param formula Object of class 'formula'; a symbolic description of
#' the model to be fitted
#' @param data Optional data frame containing the variables in the model.
#' @param waves An integer vector which identifies components in clusters.
#' The length of \code{waves} should be the same as the number of
#' observations. \code{waves} is automatically generated if none is supplied,
#' but when using \code{subset} parameter, the \code{waves} parameter must be
#' provided by the user for proper calculation.
#' @param contrasts An optional list provided when using a formula.
#' similar to \code{contrasts} from \code{glm}.
#' See the \code{contrasts.arg} of \code{model.matrix.default}.
#' @param subset An optional vector specifying a subset of observations to be
#' used in the fitting process.
#' 
#' @export
#' @rdname bisee
bisee.formula <- function(formula, data=list(),
                          clusterID,
                          waves = NULL,
                          lambda1, lambda2 = 1-lambda1,
                          contrasts = NULL,
                          subset,
                          ...)
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "clusterID", "waves", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    if (is.null(mf$clusterID)){
        mf$clusterID <- as.name("clusterID")
    }

    if (is.null(mf$waves)){
        mf$waves <- NULL
    }
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    clusterID <- model.extract(mf, clusterID)
    if(is.null(clusterID)){
        stop("clusterID variable not found.")
    }
    waves <- model.extract(mf, waves)
    
    y <- model.response(mf, "numeric")
    x <- model.matrix(attr(mf, "terms"), data=mf, contrasts)    

    ## bisee deterimes intercept based on 'intercept' parameter
    if(all(x[,1] ==1)){
        x <- x[,-1]
    }

    if(any(colSums(x) == 0)){

        cat("######## ERROR! ########\n")
        cat(colnames(x)[colSums(x) == 0])
        cat("\n")
        stop("The above factors are not found in the given observations")
    }

    
    results <- bisee.default(y, x,
                             clusterID = clusterID,
                             waves = waves,
                             lambda1 = lambda1,
                             lambda2 = lambda2,
                             ...)
    results$call <- match.call()
    results$contrasts <- contrasts
    results
}

#' @export
#' @rdname bisee
bisee.default <- function(y, x,
                          waves = NULL,
                          lambda1, lambda2 = 1-lambda1,
                          ...){



    if( lambda2 == 0){
        results <- gsee(y,x, waves = waves, ...)
    } else if(lambda1 == 0){
        results <- see(y,x, waves = waves, ...)
    }else {
        results <- bisee.fit(y, x,
                             waves = waves,
                             lambda1 = lambda1,
                             lambda2 = lambda2, ...)
    }
    results$call <- match.call()
    results
}


#' @export
#' @rdname bisee 
bisee.fit <-
function(y, x, family,
         clusterID, waves = NULL, groupID,
         corstr="independence", alpha = NULL,
         lambda1 = .5, lambda2 = 1-lambda1,
         intercept = TRUE,
         offset = 0,
         control = sgee.control(maxIt = 200, epsilon = 0.05, 
                    stoppingThreshold =  min(length(y), ncol(x))-intercept,
                                       undoThreshold = 0.005),
         standardize = TRUE,
         verbose = FALSE,
         ...){
    #######################
    ## Preliminaries/set up
    #######################

    maxIt <- control$maxIt
    epsilon <- control$epsilon
    undoThreshold <- control$undoThreshold
    interceptLimit <- control$interceptLimit
    ##If the undoThreshold is >= epsilon
    ## the check wil always trigger;
    ## so this check is added to prevent
    ## an infinte loop
    if(undoThreshold >= epsilon){
        if(verbose){
            cat(paste0("****** undoThreshold too large! reducing threshold now **********\n"))
        }
        undoThreshold  <- epsilon/10
    }
    
    if(is.null(control$stoppingThreshold)){
        stoppingThreshold <- min(length(y), ncol(x))-intercept
    } else {
        stoppingThreshold <- control$stoppingThreshold
    }
    
    p <- length(x[1,])

    if (is.character(family)){ 
        family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)){ 
        family <- family()
    }
    if(standardize){
        unstandardizedX <- x
        if(intercept){
            x <- scale(x)
        } else{
            x <- scale(x, center = FALSE)
        }
    }

    if(is.null(waves)){
        clusz <- unlist(sapply(unique(clusterID), function(x) {sum(clusterID == x)}))
        waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
    }

    
    ## currently assuming only intercept in estimation of
    ## of correlation and dispersion
    r <- 1 # number of covariates in dispersion modeling
    q <- 1 # number of covariates in correlation modeling

    ## Initial estimates for all parameters
    beta <- rep(0,p)
    phi <- stats::sd(y)^2   
    ## current initial estimation for correlation parameter
    if(is.null(alpha)){
        alpha <- 0
    }
    ## Intercept value
    beta0 <- 0

    
    meanLink <- family$linkfun
    meanLinkInv <- family$linkinv
    varianceLink <- family$variance
    mu.eta <- family$mu.eta
    
    ## Paths of parameter estimates
    path <- matrix(rep(0,(maxIt)*(p + intercept)), nrow = maxIt)
    ## currently assuming only intercept in estimating dispersion
    phiPath <- matrix(rep(0,(maxIt)*r), nrow = maxIt)
    ## currently assuming only intercept in estimating correlation
    alphaPath <- matrix(rep(0,(maxIt)*q), nrow = maxIt)


    
    numClusters <- max(clusterID)
    maxClusterSize <- max(table(clusterID))
    mu <- rep(0,length(y))

    ## it is assumed that the penalty is applied only to a single group
    ## not all ofthe covariates at once
    penalty <- function(delta,
                        lambda1,
                        lambda2,
                        scale = 1,
                        epsilon = 0){
      
      return(sqrt(length(delta))*lambda1 * sqrt(sum((scale*delta)^2)) + lambda2 * sum(abs(scale*delta)) - epsilon)
    }


    ##stoppedOn added to keep track of when the algorithm stops
    ## it assumes it goes the whole lenght unless stopped prematurely
    stoppedOn <- maxIt

    # Working correlation matrix
    R <- genCorMat(corstr = corstr, rho = alpha, maxClusterSize = maxClusterSize)
    RInv <- solve(R)
    ##################
    ## Main Algorithim
    ##################
    cat("\n")
    oldDelta <- rep(0, length(beta)) 
    it <- 0
    while (it <maxIt){
        it <- it +1
        if(verbose){
            cat(paste0("****** Beginning iteration # ", it, " **********\n"))
        }
        GEEValues <- evaluateGEE(y = y,
                                 x = x,
                                 beta = beta,
                                 beta0 = beta0,
                                 intercept,
                                 phi = phi,
                                 offset = offset,
                                 RInv = RInv,
                                 numClusters = numClusters,
                                 clusterID = clusterID,
                                 waves = waves,
                                 meanLinkInv = meanLinkInv,
                                 mu.eta = mu.eta,
                                 varianceLink = varianceLink,
                                 corstr = corstr,
                                 maxClusterSize = maxClusterSize,
                                 interceptLimit = interceptLimit)

        ## Update Estimates
        beta0 <- GEEValues$beta0
        phi <- GEEValues$phiHat
        alpha <- GEEValues$rhoHat
        RInv <- GEEValues$RInv
        
        ## Current Values of estimating Equations
        sumMean <- GEEValues$sumMean
        
        
        ###################
        ## Update Selection
        ###################
        a <- rep(0, length(unique(groupID)))
       
        for (j in unique(groupID)){
            currentGroupIndex <- groupID == j
            aCurrent <- sumMean[currentGroupIndex]

            ## extendInt may cause issues for older versions of R
            rootValues <- stats::uniroot(deltaFinder, a = aCurrent ,
                                         lambda1 =lambda1, lambda2 = lambda2,
                                         lower = 1e-100, upper = 500,
                                         extendInt = "yes")


            
            a[j] <- rootValues$root
        }
        
        ## Identify optimal group
        delta <- which(a == max(a))
        if(length(delta) > 1){
            delta <- delta[1]
        }
        theGroup <- groupID == delta
        
        ## find gamma that yields best delta for the given group
        result <- a[delta]
        groupDelta <- deltaValue(sumMean[theGroup], lambda1, lambda2, result)

        ## Scaling factor to make step of appropriate size
        c <- epsilon/(lambda2*sum(abs(groupDelta)) + lambda1*sqrt(length(groupDelta)))
        ## final answer for a given group
        groupDelta <- groupDelta * c

        ## Check if the update is effectively undoing the last one
        if (sum(abs(oldDelta[theGroup] + groupDelta))<= undoThreshold){
            if(verbose){
                cat(paste0("****** Step Undone! Reducing Stepsize **********\n"))
            }

            if(it>2){
                if (intercept){
                    beta <- path[it - 2,-1]
                } else {
                    beta <- path[it - 2,]
                }
            } else{
                beta <- rep(0, length(beta))
            }
            
            epsilon <- epsilon/2
            it <- it - 2
            oldDelta <- rep(0, length(beta))

            ##If the undoThreshold is >= epsilon
            ## the check wil always trigger;
            ## so this check is added to prevent
            ## an infinte loop
            if(undoThreshold >= epsilon){
                if(verbose){
                    cat(paste0("****** undoThreshold too large! reducing threshold now **********\n"))
                }
                undoThreshold  <- epsilon/10
            }
            
            ## If the check is passed and the update
            ## is sufficiently different from the previous
            ## update
        } else {
            oldDelta <- rep(0, length(beta))
            oldDelta[theGroup] <- groupDelta


            ## Update estimates
            beta[theGroup] <- beta[theGroup] +  groupDelta

            ## Update paths
            if(intercept){
                path[it,] <- c(beta0, beta)
            }
            else{
                path[it,] <-  beta
            }
            phiPath[it,] <- phi
            alphaPath[it,] <- alpha

            ###########
            ## stopping mechanism when the alogrithim has reached saturation
            if((sum(beta != 0) >= stoppingThreshold) & (it< maxIt) ){
                print("stopped on")
                print(it)
                print(a[delta])
                path[((it+1):maxIt),] <- matrix(rep(path[it,],(maxIt - it)),
                                                nrow = (maxIt - it),
                                                byrow = TRUE)
                phiPath[((it+1):maxIt),] <- matrix(rep(phiPath[it,],(maxIt - it)),
                                                   nrow = (maxIt - it),
                                                   byrow = TRUE)
                alphaPath[((it+1):maxIt),] <- matrix(rep(alphaPath[it,],(maxIt - it)),
                                                     nrow = (maxIt - it),
                                                     byrow = TRUE)
                
                stoppedOn <- it
                break
            }

        }
    }
    
    result <- list(path = path,
                   gammaPath = phiPath,
                   alphaPath = alphaPath,
                   stoppedOn = stoppedOn,
                   maxIt = maxIt,
                   y = y,
                   x = x,
                   intercept = intercept,
                   clusterID = clusterID,
                   groupID = groupID,
                   family = family,
                   offset = offset,
                   epsilon = epsilon)

    if(standardize){
        result$x <- unstandardizedX
        temp <- path
        
        if(intercept){
            temp[,-1] <- t(t(path[,-1]) /attr(x, "scaled:scale"))
            temp[,1] <- path[,1] - crossprod(t(path[,-1]),
                                             (attr(x, "scaled:center")/attr(x, "scaled:scale")))  
        } else{
            temp <- t(t(path) /attr(x, "scaled:scale"))
        }
        
        result$standardizedPath <-  path
        result$path <- temp
        result$standardizedX <- x
    }

    
    class(result) <- "sgee"

    result
}
