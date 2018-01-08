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
#################################################################################' Stagewise Estimating Equations Implementation
#'
#' 
#' Function to perform SEE, a  Forward Stagewise regression
#' approach for model selection / dimension reduction
#' using Generalized Estimating Equations
#' 
#' Function to implement SEE, a stagewise regression approach
#' that is designed to perform model selection in the context of
#' Generalized Estimating Equations. Given a response \code{y} and
#' a design matrix \code{x}
#' (excluding intercept) SEE generates a path of stagewise regression
#' estimates for each covariate based on the provided step size \code{epsilon}.
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
#' A stochastic version of this function can also be called. using the
#' auxiliary function \code{sgee.control} the parameters \code{stochastic},
#' \code{reSample}, and \code{withReplacement} can be given to \code{see}
#' to perform a sub sampling step in the procedure to make the SEE
#' implementation scalable for large data sets.
#' 
#'
#' @param y Vector of response measures that corresponds with modeling family
#' given in 'family' parameter. \code{y} is assumed to be the same length as
#' \code{clusterID} and is assumed to be organized into clusters as dictated by
#' \code{clusterID}.
#' @param x Design matrix of dimension \code{length(y)} x nvar,
#' the number of variables, where each row is
#' represents an observation of predictor variables. 
#' @param family Modeling family that describes the marginal distribution of
#' the response. Assumed to be an object such as \code{gaussian()} or
#' \code{poisson()}.
#' @param clusterID Vector of integers that identifies the clusters of response
#' measures in \code{y}. Data and \code{clusterID} are assumed to
#' 1) be of equal lengths,
#' 2) sorted so that observations of a cluster are in contiguous rows, and 3)
#' organized so that \code{clusterID} is a vector of consecutive integers.
#' @param corstr A character string indicating the desired working correlation
#' structure. The following are implemented : "independence" (default value),
#' "exchangeable", and "ar1".
#' @param alpha An initial guess for the correlation parameter value
#' between -1 and 1 . If left NULL (the default), the initial estimate is 0.
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
#' estimates, the iteration at which SEE terminated, and initial regression
#' values including \code{x}, \code{y}, code{family}, \code{clusterID},
#' \code{groupID}, \code{offset}, \code{epsilon}, and \code{numIt}.
#' 
#' @author Gregory Vaughan
#' @references Vaughan, G., Aseltine, R., Chen, K., Yan, J., (2017). Stagewise
#' Generalized Estimating Equations with Grouped Variables. Biometrics 73,
#' 1332-1342. URL: http://dx.doi.org/10.1111/biom.12669,
#' doi:10.1111/biom.12669.
#'
#' Wolfson, J. (2011). EEBoost: A general method for prediction
#' and variable selection based on estimating equations. Journal of the
#' American Statistical Association 106, 296--305.
#'
#' Tibshirani, R. J. (2015). A general framework for fast stagewise
#' algorithms. Journal of Machine Learning Research 16, 2543--2588.
#'
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
#' groupSize <- 1
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
#'
#' ## Perform Fitting by providing formula and data
#' genDF <- data.frame(generatedData$y, generatedData$x)
#' names(genDF) <- c("Y", paste0("Cov", 1:p))
#' coefMat1 <- see(formula(genDF), data = genDF,
#'                  family = gaussian(),
#'                  waves = rep(1:4, 50), 
#'                  clusterID = generatedData$clusterID,
#'                  groupID = generatedData$groupID, 
#'                  corstr = "exchangeable",
#'                  control = sgee.control(maxIt = 50, epsilon = 0.5),
#'                  verbose = TRUE)
#'
#' ## set parameter 'stochastic' to 0.5 to implement the stochastic
#' ## stagewise approach where a subsmaple of 50% of the data is taken
#' ## before the path is calculation.
#' ## See sgee.control for more details about the parameters for the
#' ## stochastic stagewise approach
#' 
#' coefMat2 <- see(formula(genDF), data = genDF,
#'                  family = gaussian(),
#'                  waves = rep(1:4, 50), 
#'                  clusterID = generatedData$clusterID,
#'                  groupID = generatedData$groupID, 
#'                  corstr = "exchangeable",
#'                  control = sgee.control(maxIt = 50, epsilon = 0.5,
#'                                         stochastic = 0.5), 
#'                  verbose = FALSE)
#' 
#' par(mfrow = c(2,1))
#' plot(coefMat1)
#' plot(coefMat2)
#' 
#' @export see
#' @name see
NULL

#' @export
#' @rdname see 
see <- function(y, ...) UseMethod("see")

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
#' @rdname see
see.formula <- function(formula, data=list(),
                          clusterID,
                          waves = NULL,
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

    ## see determines intercept based on 'intercept' parameter
    if(all(x[,1] ==1)){
        x <- x[,-1]
    }

    if(any(colSums(x) == 0)){

        cat("######## ERROR! ########\n")
        cat(colnames(x)[colSums(x) == 0])
        cat("\n")
        stop("The above factors are not found in the given observations")
    }

    
    results <- see.default(y, x,
                           clusterID = clusterID,
                           waves = waves,
                           ...)
    results$call <- match.call()
    results$contrasts <- contrasts
    results
}

#' @export
#' @rdname see
see.default <- function(y, x,
                          waves = NULL,
                          ...){




    results <- see.fit(y,x, waves = waves, ...)
    results$call <- match.call()
    results
}



#' @export
#' @rdname see
see.fit <- function(y, x,  family,
                    clusterID,
                    waves = NULL,
                    corstr="independence", alpha = NULL,
                    intercept = TRUE,
                    offset = 0,
                    control = sgee.control(maxIt = 200, epsilon = 0.05, 
                        stoppingThreshold =  min(length(y), ncol(x))-intercept,
                        undoThreshold = 0),
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

    ##############################
    ## Stochastic parameter values
    ##############################
    ## stochastic == 1 defaults to normal stagewise approach
    ## with no subsampling
    stochastic = control$stochastic 
    sampleProb = control$sampleProb
    reSample = control$reSample
    withReplacement = control$withReplacement
    
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
    
    p <- ncol(x)

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
        ## minor alteration that allows
        ## for columns that don't vary at all
        x[,attr(x, "scaled:scale") ==0] <- 0
        attr(x, "scaled:scale")[attr(x, "scaled:scale") ==0] <- 1
    }

    if(is.null(waves)){
        clusz <- unlist(sapply(unique(clusterID), function(x) {sum(clusterID == x)}))
        waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
    }
    
    ## currently assuming only intercept in estimation
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
    colnames(path) <- c("(Intercept)",
                        colnames(unstandardizedX))[c(intercept, rep(TRUE, p))]
    ## currently assuming only intercept in estimating dispersion
    phiPath <- matrix(rep(0,(maxIt)*r), nrow = maxIt)
    ## currently assuming only intercept in estimating correlation
    alphaPath <- matrix(rep(0,(maxIt)*q), nrow = maxIt)

    clusterIDs <- unique(clusterID)
    numClusters <- length(clusterIDs)
    maxClusterSize <- max(waves)

    if(stochastic<1){
        ## determine size of sub-sample
        ## ceiling is used as round may be inconsistent
        ## from OS to OS
        sampleSize <- ceiling(numClusters*stochastic)
        
        ## stochastic == 1 indicates the Deterministic approach
    } else if(stochastic == 1){
        yCurr <- y
        xCurr <- x
        clusterIDCurr <- clusterID
        wavesCurr <- waves
    }

    ##stoppedOn added to keep track of when the algorithm stops
    ## it assumes it goes the whole length unless stopped prematurely
    stoppedOn <- maxIt

    ## Working correlation matrix
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

        #########################
        ## Stochastic Subsampling
        #########################
        ## stochastic <1  --> stochastic == 1 implies deterministic
        ## approach, no subsampling
        ## it == 1 --> if we are subsampling, we always subsample
        ## the first iteration
        ## (reSample > 0 & ((it+1) %% reSample == 0)) --> if
        ## reSample >0, then that means we are sampling at
        ## regular intervals; reSample == 1 implies a subsample
        ## is collected very step, reSample == 2 implies a subsample
        ## is collected every other steps, etc ...
        if(stochastic <1 & (it == 1 | (reSample > 0 & ((it+1) %% reSample == 0)))){
            ## Calculate the actual subsampling distirbution.
            ## if a function was given to determine the subsampling
            ## probabilities, then that distribution must be
            ## calculated every time a sub-sample is collected
            if(is.function(sampleProb) | it == 1){                
                ## function to initialize sampling Distribtuion
                sampleDist <- samplingDistCalculation(sampleProb = sampleProb,
                                                      y = y,
                                                      x = x,
                                                      clusterID = clusterID,
                                                      waves = waves,
                                                      beta = beta,
                                                      beta0 = beta0,
                                                      phi = phi,
                                                      alpha = alpha,
                                                      offset = offset,
                                                      meanLinkInv = meanLinkInv,
                                                      varianceLink = varianceLink,
                                                      corstr = corstr,
                                                      mu.eta = mu.eta)
                
            }
            
            subsampleOutput <- subsample(sampleDist = sampleDist,
                                         sampleSize = sampleSize,
                                         withReplacement = withReplacement,
                                         clusterIDs = clusterIDs,
                                         clusterID = clusterID)
            

            subSampleIndicator <- subsampleOutput$subSampleIndicator
            ## if we sample with replacement, the clusterIDCurr
            ## variable has to be constructed differently than if we
            ## sample without replacement
            clusterIDCurr <- subsampleOutput$clusterIDCurr
            
            yCurr <- y[subSampleIndicator]
            xCurr <- x[subSampleIndicator,]
            wavesCurr <- waves[subSampleIndicator]
            
            ## if the sub-sample proportion is too small,
            ## it is possible
            ## that in a given iteration there are not
            ## enough observations
            ## for the number of features in the model.
            ## To prevent this case,
            ## the stopping threshold is adjusted here
            ## if it is higher than the number of observations
            stoppingThreshold <- min(stoppingThreshold, length(yCurr)-intercept)
            
            
        } 
        
         GEEValues <- evaluateGEE(y = yCurr,
                                 x = xCurr,
                                 beta = beta,
                                 beta0 = beta0,
                                 intercept,
                                 phi = phi,
                                 offset = offset,
                                 RInv = RInv,
                                 numClusters = sampleSize,
                                 clusterID = clusterIDCurr,
                                 waves = wavesCurr,
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

        ## Current values of the estimating equations
        sumMean <- GEEValues$sumMean

        ## Identify optimal update
        a <- abs(sumMean)   
        delta <- which(a== max(a))

        if(verbose){
            print(paste0("L1 Norm of Gradient :", sum(a)))
            print(paste0("updating index: ", delta, " by ", epsilon))
        }

        ## Check if the update is effectively undoing the last one
        if (sum(abs(oldDelta[delta] + epsilon*sign(sumMean[delta])))<= undoThreshold){
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
            oldDelta[delta] <- epsilon*sign(sumMean[delta])

            ## Update estimates
            beta[delta] <- beta[delta] +epsilon*sign(sumMean[delta]) 

            ## Update the paths
            if(intercept){
                path[it,] <- c(beta0, beta)
            }
            else{
                path[it,] <-  beta
            }
            phiPath[it,] <- phi
            alphaPath[it,] <- alpha

            if(verbose){
                print("current and previous steps")
                subPath <- path
                dimnames(subPath) <- list(1:nrow(path), 1:ncol(path))
                print(subPath[c(it-1, it), colSums(path) !=0])
            }
        
            ###########
            ## stopping mechanism when the alogrithim has
            ## reached saturation.
            ## sum(a) <0.5 threshold added to prevent possible loop
            ## that can happen with binary data using the adaptive
            ## step size where the algorithm thinks a step keeps
            ## being undone, but really the estimating equations
            ## are all VERY close to 0
            if(((sum(beta != 0) >= stoppingThreshold) | sum(a) < 0.5 )& (it< maxIt) ){
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
                   groupID = 1:ncol(x),
                   family = family,
                   offset = offset,
                   epsilon = epsilon,
                   stochastic = stochastic, 
                   sampleProb = sampleProb,
                   reSample = reSample,
                   withReplacement = withReplacement)

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
