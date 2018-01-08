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
#' Interaction stagewise estimating equations
#' 
#' Perform model selection with clustered data while considering interaction
#' terms using one of two stagewise methods. The first (ACTS) uses an active set
#' approach in which interaction terms are only considered for a given update
#' if the corresponding main effects have already been added to the model.
#' The second approach (HiLa) approximates the regularized path for
#' hierarchical lasso with Generalized Estimating Equations. In this second
#' approach, the model hierarchy is guaranteed in each individual step, thus
#' ensuring the desired hierarchy throughout the path.
#' 
#' @param y Vector of response measures that corresponds with modeling family
#' given in 'family' parameter. 'y' is assumed to be the same length as
#' 'clusterID' and is assumed to be organized into clusters as dictated by
#' 'clusterID'.
#' @param x Design matrix of dimension length(y) x nvars where each row is
#' represents an obersvation of predictor variables. Assumed to be scaled.
#' @param family Modeling family that describes the marginal distribution of
#' the response. Assumed to be an object such as 'gaussian()' or 'poisson()'
#' @param clusterID Vector of integers that identifies the clusters of response
#' measures in 'y'. Data and 'clusterID' are assumed to 1) be of equal lengths,
#' 2) sorted so that observations of a cluster are in contiguous rows, and 3)
#' organized so that 'clusterID' is a vector of consecutive integers.
#' @param interactionID A (p^2+p)/2 x 2 matrix of interaction IDs. Main effects
#' have the same (unique) number in both columns for their corresponding row.
#' Interaction effects have each of their corresponding main effects in the
#' two columns. it is assumed that main effects are listed first. It is
#' assumed that the main effect IDs used start at 1 and go up tp the number
#' of main effects, p.
#' @param corstr A character string indicating the desired working correlation
#' structure. The following are implemented : "independence" (default value),
#' "exchangeable", and "ar1".
#' @param alpha An intial guess for the correlation parameter value
#' between -1 and 1 . If left NULL (the default), the initial estimate is 0.
#' @param intercept Binary value indicating where an intercept term is
#' to be included in the model for estimation. Default is to include an
#' intercept.
#' @param offset Vector of offset value(s) for the linear predictor. 'offset'
#' is assumed to be either of length one, or of the same length as 'y'.
#' Default is to have no offset.
#' @param method A character string indicating desired method to be used to
#' perform interaction selection. Value can either be "ACTS", where an active
#' set approach is taken and interaction terms are considered for selection
#' only after main effects are brought in, or "HiLa", where the hierarchical
#' lasso penalty is used to ensure hierarchy is maintained in each step.
#' Default Value is "ACTS".
#' @param control A list of parameters used to contorl the path generation
#' process; see \code{sgee.control}.
#' @param standardize A logical parameter that indicates whether or not
#' the covariates need to be standardized before fitting (but after generating
#' interaction terms from main covariates).
#' If standardized before fitting, the unstandardized
#' path is returned as the default, with a \code{standardizedPath} and
#' \code{standardizedX} included
#' separately. Default value is \code{TRUE}.
#' @param verbose Logical parameter indicating whether output should be produced
#' while isee is running. Default value is FALSE.
#' @param ... Not currently used
#' 
#' @return Object of class 'sgee' containing the path of coefficient estimates,
#' the path of scale estimates, the path of correlation parameter
#' estimates, and the iteration at which iSEE terminated, and initial regression
#' values including \code{x}, \code{y}, code{family}, \code{clusterID},
#' \code{interactionID}, \code{offset}, \code{epsilon}, and \code{numIt}.
#' 
#' @note While the two different possible methods that can be used with
#' \code{isee} reflect two different "styles" of stagewise estimation,
#' both achieve a desired hierarchy in the resulting model paths.
#'
#' When considering models with interaction terms, there are three forms
#' of hierarchy that may be present. Strong hierarchy implies that
#' interaction effects are included in the model only if both of its
#' corresponding main effects are also included in the model. Weak hierarchy
#' implies that an interaction effect can be in the model only if AT LEAST
#' one of its corresponding main effects is also included. The third type
#' of hierarchy is simply a lack of hierarchy; that is an interaction term
#' can be included regardless of main effects.
#'
#' In practice strong hierarchy is usually what is desired as it is the
#' simplest to interpret, but requires a higher amount of computation when
#' performing model selection. Weak hierarchy is sometimes used as a compromise
#' between the interpret-ability of strong hierarchy and the computational ease
#' of no hierarchy. Both \code{isee} methods only implement strong hierarchy
#' as the use of stagewise procedures greatly reduces the computational burden.
#'
#' The active set appraoch, ACTS, tends to have slightly better predictive
#' and model selection performance when the true model is closer to a purely
#' strong hierarchy, but HiLa tends to do better if the true model hierarchy
#' is closer to having a purely weak hierarchy. Thus, in practice, it is
#' important to use external information and judgement to determine which
#' approach is more appropriate.
#' 
#' @author Gregory Vaughan
#' @references Vaughan, G., Aseltine, R., Chen, K., Yan, J., (2017). Efficient
#' interaction selection for clustered data via stagewise generalized
#' estimating equations.  Department of Statistics, University of
#' Connecticut. Technical Report.
#'
#' Zhu, R., Zhao, H., and Ma, S. (2014). Identifying
#' gene-environment and gene-gene interactions using a progressive
#' penalization approach. Genetic Epidemiology 38, 353--368.
#' @references Bien, J., Taylor, J., and Tibshirani, R. (2013). A lasso
#' for hierarchical interactions. The Annals of Statistics 41, 1111--1141.
#' 
#' @examples
#' 
#' #####################
#' ## Generate test data
#' #####################
#' 
#' ## Initialize covariate values
#' p <- 5 
#' beta <- c(1, 0, 1.5, 0, .5, ## Main effects
#'           rep(0.5,4), ## Interaction terms
#'           0.5, 0, 0.5,
#'           0,1,
#'           0)
#' 
#' 
#' generatedData <- genData(numClusters = 50,
#'                          clusterSize = 4,
#'                          clusterRho = 0.6,
#'                          clusterCorstr = "exchangeable",
#'                          yVariance = 1,
#'                          xVariance = 1,
#'                          beta = beta,
#'                          numMainEffects = p,
#'                          family = gaussian(),
#'                          intercept = 1)
#' 
#'  
#' ## Perform Fitting by providing formula and data
#' genDF <- data.frame(Y = generatedData$y, X = generatedData$xMainEff)
#'
#' ## Using "ACTS" method
#' coefMat1 <- isee(formula(paste0("Y~(",
#'                                paste0("X.", 1:p, collapse = "+"),
#'                                  ")^2")),
#'                   data = genDF,
#'                   family = gaussian(),
#'                   clusterID = generatedData$clusterID,
#'                   corstr = "exchangeable",
#'                   method = "ACTS",
#'                   control = sgee.control(maxIt = 50, epsilon = 0.5))
#'
#' ## Using "HiLa" method
#' coefMat2 <- isee(formula(paste0("Y~(",
#'                                paste0("X.", 1:p, collapse = "+"),
#'                                  ")^2")),
#'                   data = genDF,
#'                   family = gaussian(),
#'                   clusterID = generatedData$clusterID,
#'                   corstr = "exchangeable",
#'                   method = "HiLa",
#'                   control = sgee.control(maxIt = 50, epsilon = 0.5))
#' 
#' @export isee
#' @name isee
NULL


#' @export 
#' @rdname isee 
isee <- function(y, ...) UseMethod("isee")

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
#' @rdname isee
#' @export
isee.formula <- function(formula, data=list(),
                         clusterID,
                         waves = NULL,
                         interactionID = NULL,
                         contrasts = NULL,
                         subset,
                         method = "ACTS",
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


    ## isee deterimes intercept based on 'intercept' parameter
    if(all(x[,1] ==1)){
        x <- x[,-1, drop = FALSE]
    }
    
    ## Making the interaction ID if not provided
    
    if(is.null(interactionID)){
        ## using mf directly to get the factorMatrix
        ## does not yield the desired matrix
        ## instead, the following workaround is used
        tempFormulaString <- paste0("y ~ ",
                                    paste0(gsub(pattern = ":",
                                                replacement = "*",
                                                colnames(x)),
                                           collapse = " + "))
        tempDF <- data.frame(cbind(y,x)) 
        tempModelFrame <- stats::model.frame(as.formula(tempFormulaString),
                                             data = tempDF)

        factorMatrix <- attr( attr(tempModelFrame, "terms"), "factors")

        ## apply should come out as a matrix
        interactionID <- t(apply(factorMatrix[-1,,drop=FALSE],
                                 MARGIN = c(2),
                                 FUN = function(x){
                                     index <- which(x != 0)
                                     if(length(index)==1){
                                         index <- rep(index,2)}
                                     index
                                 }))
    }

    
    



    if((ncol(interactionID) != 2) |  (ncol(x) != nrow(interactionID))){
        stop("interactionID provided not of correct dimensions. Should
be (p^2 + p)/2 by 2 matrix, where p is the number of effects being considered")
    }

    if(any(colSums(x) == 0)){

        cat("######## ERROR! ########\n")
        cat(colnames(x)[colSums(x) == 0])
        cat("\n")
        stop("The above factors are not found in the given observations")
    }
    
    results <- isee.default(y, x,
                            clusterID = clusterID,
                            waves = waves,
                            interactionID = interactionID,
                            method = method,
                            ...) 

    results$call <- match.call()
    
    results
}

#' @export
#' @rdname isee
isee.default <- function(y, x,
                         waves = NULL,
                         interactionID,
                         method = "ACTS",
                         ...){

    if(method == "ACTS"){
        results <- acts.fit(y, x,
                            waves = waves,
                            interactionID = interactionID,
                            method = method,
                            ...)
    } else if (method == "HiLa"){
        results <- hila.fit(y, x,
                             waves = waves,
                             interactionID = interactionID,
                             ...)
    } else{
        stop("paramter 'method' must be supplied a value of either 'ACTS', or 'HiLa'")
    }

    results$method <- method
    results$call <- match.call()
    results
}



#' @export
#' @rdname isee
acts.fit <-
    function(y, x,
             interactionID,
             family,
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


    ## interaction terms are marked with main effects
    ## so largest value in interactionID is the number
    ## of covariates
    p <- max(interactionID)

    interactionManager <- list()
    ## initalize activeSet to be all the main effects
    ## used to actually subsetthe estimating equations
    ## active set requires both columns ot be TRUE in order for
    ## an effecto be considered
    interactionManager$activeSet <- matrix(rep(interactionID[,1] == interactionID[,2],2),
                                           ncol = 2)

    interactionManager$index <- 1:nrow(interactionID)
    
    ## active effects keeps track of which main effects
    ## have been included, used to update activeSet
    activeEffects <-c()

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

    ## Initalize estimates for all parameters
    ## beta has to match X, which may have more than p columns
    beta <- rep(0,ncol(x))
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
    ## has to match number of covariates, which may be larger than p
    path <- matrix(rep(0,(maxIt)*(ncol(x) + intercept)), nrow = maxIt)
    ## currently assuming only intercept in estimating dispersion
    phiPath <- matrix(rep(0,(maxIt)*r), nrow = maxIt)
    ## currently assuming only intercept in estimating correlation
    alphaPath <- matrix(rep(0,(maxIt)*q), nrow = maxIt)


    clusterIDs <- unique(clusterID)
    numClusters <- length(clusterIDs)
    maxClusterSize <- max(waves)

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



        ## interactionManager$activeSet[,1] & interactionManager$activeSet[,2]
        ## is meant to keep track of what effects are part of
        ## the active set

        ## subset estimating equations, only look at
        ## activeSet
        a <- abs(sumMean * (interactionManager$activeSet[,1] & interactionManager$activeSet[,2]))
        
        ## Identify optimal group
        delta <- which(a == max(a))

        aCurrent <- sumMean[delta]
        
        
        
        ## Check if the update is effectively undoing the last one
        if (sum(abs(oldDelta[delta] + epsilon * sign(aCurrent)))<= undoThreshold){
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
            oldDelta[delta] <- epsilon * sign(aCurrent)


            ## Update estimate
            beta[delta] <- beta[delta] +  epsilon * sign(aCurrent)

 

            ## Update activeSet
            ## first check if a main effect was added
            theEffect <- interactionManager$index[delta]
            if(theEffect <= p){
                ## then check if the activeEffects needs to be updated
                if(!(theEffect %in% activeEffects)){
                    ## update the active effects
                    activeEffects <- c(activeEffects, theEffect)
                    interactionManager$activeSet <- interactionManager$activeSet  | interactionID == theEffect
                }
            }

        
            ## Update Paths
            if(intercept){
                path[it,] <- c(beta0, beta)
            }
            else{
                path[it,] <-  beta
            }
            phiPath[it,] <- phi
            alphaPath[it,] <- alpha
            
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
                   interactionID = interactionID,
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


#' @export
#' @rdname isee
hila.fit <-
    function(y, x,
             interactionID,
             family,
             clusterID,
             waves = NULL,
             corstr="independence", alpha = NULL, 
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
    

    ## interaction terms are marked with main effects
    ## so largest value in interactionID is the number
    ## of covariates
    p <- max(interactionID)

    interactionManager <- list()
    ## initalize activeSet to be all the main effects
    ## used to actually subsetthe estimating equations
    ## active set requires both columns ot be TRUE in order for
    ## an effecto be considered
    interactionManager$activeSet <- matrix(rep(interactionID[,1] == interactionID[,2],2),
                                           ncol = 2)
    interactionManager$mainEffects <- interactionManager$activeSet[,1]
    interactionManager$index <- 1:nrow(interactionID)
    
    ## active effects keeps track of which main effects
    ## have been included, used to update activeSet
    activeEffects <-c()

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

    ## Initalize estimates for all parameters
    ## beta has to match X, which may have more than p columns
    beta <- rep(0,ncol(x))
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
    ## has to match number of covariates, which may be larger than p
    path <- matrix(rep(0,(maxIt)*(ncol(x) + intercept)), nrow = maxIt)
    ## currently assuming only intercept in estimating dispersion
    phiPath <- matrix(rep(0,(maxIt)*r), nrow = maxIt)
    ## currently assuming only intercept in estimating correlation
    alphaPath <- matrix(rep(0,(maxIt)*q), nrow = maxIt)


    clusterIDs <- unique(clusterID)
    numClusters <- length(clusterIDs)
    maxClusterSize <- max(waves)

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

        ## subset estimating equations, only look at
        ## activeSet

        ## for now disabling the active set funcitonality
        ##currentActive <- interactionManager$activeSet[,1] & interactionManager$activeSet[,2]
        ##currentIndecies <- interactionManager$index[currentActive]
        ##a <- abs(sumMean[currentActive])
        a <- abs(sumMean)
        

        ## aFull is the vector of 3|U_{ii}| (for main effects)
        ## and |U_{ij}| +|U_{ii}|+|U_{jj}| (for interaction effects)
        aFull <- rowSums(matrix(c(a, c(a)[interactionID]) ,ncol = 3))
        fullMax <- max(aFull)
        delta <- (aFull == fullMax)

        effectIndecies <- interactionID[delta,]
        if(effectIndecies[1] != effectIndecies[2]){
            effect1 <- interactionID == effectIndecies[1]
            effect2 <- interactionID == effectIndecies[2]

            theEffect <- (rowSums(effect1) == 2) | delta | (rowSums(effect2) == 2)
        } else{
            theEffect <- delta
        }
        

        aCurrent <- sumMean[theEffect]
        
                ## Check if the update is effectively undoing the last one
        if (sum(abs(oldDelta[theEffect] + epsilon * sign(aCurrent)/ sum(theEffect)))<= undoThreshold){
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
            oldDelta[theEffect] <- epsilon * sign(aCurrent)/ sum(theEffect)
            
        
            ## Update estimate
            ## sum(theEffect) is either 1 or 3
            beta[theEffect] <- beta[theEffect] +  epsilon * sign(aCurrent)/ sum(theEffect)
            ## alternate update form not validated in paper
            ## beta[theEffect] <- beta[theEffect] +  epsilon * aCurrent/ sum(abs(aCurrent))

            ## Update activeSet
            ## first check if a main effect was added
        
            ## for now disabling the active set funcitonality
            ## if(theEffect <= p){
            ## then check if the activeEffects needs to be updated
            ##    if(!(theEffect %in% activeEffects)){
            ## update the active effects
            ##        activeEffects <- c(activeEffects, theEffect)
            ##        interactionManager$activeSet <- interactionManager$activeSet  | interactionID == theEffect
            ##    }
            ##}

            
            ## Update Paths
            if(intercept){
                path[it,] <- c(beta0, beta)
            }
            else{
                path[it,] <-  beta
            }
            phiPath[it,] <- phi
            alphaPath[it,] <- alpha

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
                   interactionID = interactionID,
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
    
