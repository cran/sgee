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
#################################################################################' Response and  Covariate Data Generation
#' 
#' Function to generate data that can be used to test Forward stagewise /
#' Penalized Regression techniques. Currently marginally Gaussian and
#' Poisson responses are possible.
#'  
#' Function is provided to allow the user simple data generation as
#' \code{sgee} functions were designed for.
#' Various parameters controlling
#' aspects such as the response correlation, the covariate group
#' structure, the marginal response distribution, and the signal to
#' noise ratio for marginally gaussian responses are
#' provided to allow a great deal of specificity over the kind of data
#' that is generated. 
#' 
#' @param numClusters Number of clusters to be generated.
#' @param clusterSize Size of each cluster.
#' @param clusterRho Correlation parameter for response.
#' @param clusterCorstr String indicating cluster Correlation structure.
#' Parameter is fed to \code{genCorMat},
#' so all possible entries for \code{genCorMat} are allowed.
#' @param yVariance Optional scalar value specifying the marginal response
#' variance; overrides \code{SNR}. 
#' @param xVariance Scalar value indicating marginal variance of the
#' covariates.
#' @param numGroups Number of covariate groups to be generated. Default
#' behavior is to generate groups of size 1 (effectively no groups).
#' If covariate groups are desired, \code{numGroups} and \code{groupSize}
#' must be given such that \code{length(beta)} equals
#' \code{numGroups} * \code{groupSize}.
#' @param groupSize Size of each group.
#' @param groupRho Within group correlation parameter.
#' @param beta Vector of coefficient values used to generate response.
#' @param numMainEffects An integer indicating that the first
#' \code{numMainEffects} terms in \code{beta} are to be treated
#' as main effects and the remaining terms are pairwise interaction effects,
#' which are in the same order as generated by \code{model.matrix}.
#' Default value of \code{NULL} indicates no interaction terms are included.
#' The use of \code{numMainEffects} overrides any covariate grouping structure
#' provided by the user.
#' @param family Marginal response family; currently \code{gaussian()} and
#' \code{poisson()} are accepted.
#' @param SNR Scalar value that allows fixing the signal
#' to noise ratio as defined as the ratio of the (observed) variance in the
#' linear predictor to the variance of the response conditioned on the
#' covariates.
#' @param intercept Scalar value indicating the true intercept value.
#' @return List containing the generated response, \code{y}, the generated
#' covariates, \code{x}, a vector identifying the responses clusters,
#' \code{clusterID}, and a vector identifying the covariate groups,
#' \code{groupID}.
#' 
#' @note Function is ued to generate both the desired covariate structure and
#' the desired response structure. To generate poisson responses, functions
#' from the R package \code{coupla} are used.
#'
#' Current implementation of interactions overwrites any previous grouping
#' structure; that is the number of groups becomes p and the group sizes
#' are set to 1.
#' @author Gregory Vaughan
#' @examples
#' 
#' 
#' ## A resonse variance can be given,
#' dat1 <- genData(numClusters = 10,
#'                 clusterSize = 4,
#'                 clusterRho = .5,
#'                 clusterCorstr = "exchangeable",
#'                 yVariance = 1,
#'                 xVariance = 1,
#'                 numGroups = 5,
#'                 groupSize = 4,
#'                 groupRho = .5,
#'                 beta = c(rep(1,8), rep(0,12)),
#'                 family = gaussian(),
#'                 intercept = 1)
#' 
#' ## or the signal to noise ratio can be fixed
#' dat2 <- genData(numClusters = 10,
#'                 clusterSize = 4,
#'                 clusterRho = .5,
#'                 clusterCorstr = "exchangeable",
#'                 xVariance = 1,
#'                 numGroups = 5,
#'                 groupSize = 4,
#'                 groupRho = .5,
#'                 beta = c(rep(1,8), rep(0,12)),
#'                 family = poisson(),
#'                 SNR = 10,
#'                 intercept = 1)
#' 
#' @export genData
#' @import stats
genData <-function(numClusters,
                   clusterSize = 1,
                   clusterRho = 0,
                   clusterCorstr = "exchangeable",
                   yVariance = NULL,
                   xVariance = 1,
                   numGroups = length(beta),
                   groupSize = 1,
                   groupRho = 0,
                   beta = 0,
                   numMainEffects = NULL,
                   family = gaussian(),
                   SNR = NULL,
                   intercept = 0){

  
    n <- numClusters * clusterSize

    ZVariance <- groupRho/(1-groupRho) 

    ## rho and corstr can be changed to
    ## induce correlation between groups
    corstr <- "independence"
    rho <- 0

    if(is.null(numMainEffects)){
        p <- length(beta)
    } else{
        p <- numMainEffects
        ## for now, the use of interaction terms will override any
        ## grouping structure
        numGroups <- p
        groupSize <- 1
        mainEffID <- matrix(rep(1:p,each =2), nrow = 2)
        interactionID <- t(cbind(mainEffID, utils::combn(p,2)))
        if(length(beta) != (p^2 +p) /2){
            stop("length of beta must match number of main effects plus all pairwise interactions")
        }
    }

    if(p != groupSize * numGroups){
        stop("length of beta must equal groupSize * numGroups")
    }
    
    r <- matrix(rnorm(p* n), ncol = p)
    Z <- mvtnorm::rmvnorm(n,
                          mean = rep(0, numGroups),
                          sigma = genCorMat(corstr,
                                          rho,
                                          numGroups) * ZVariance)
    
    D <- matrix(rep(t(Z), times = rep(groupSize, numGroups*n)), nrow = n, byrow = TRUE)

    ## Induced within group exchangeable correlation
    ## with correlation of groupRho
    ## rho and corstr can be changed to induce
    ## between group correlation as well
    X <- (D + r)/sqrt(1 + ZVariance) * xVariance

    if(!is.null(numMainEffects)){
        ## Make the interaction columns
        XmainEff <- X
        Xdf <- data.frame(X)
        form <- formula(paste0("~(", paste0("X", 1:ncol(X), collapse = "+"), ")^2-1"))
        X <- model.matrix(form, Xdf)
        
    }
    
    mn.Y <- X %*% beta + intercept

    if(is.null(yVariance)){
        if(is.null(SNR)){
            stop("yVariance or signal to noise ratio must be provided")
        } else{
            yVariance <- c(var(mn.Y))/SNR
        }
    }

    clusterCorrelationMatrix <- genCorMat(clusterCorstr,
                                        clusterRho,
                                        clusterSize)

    if(family$family == "gaussian"){
        Y <- rep(NA , numClusters * clusterSize)
        for(clustNum in 1:numClusters){
            Y[((clustNum-1)*clusterSize +1):(clustNum*clusterSize)] <- mvtnorm::rmvnorm(1,
                                                                                        mean=mn.Y[((clustNum-1)*clusterSize +1):(clustNum*clusterSize)],
                                                                                        sigma= yVariance * clusterCorrelationMatrix)
        }
    } else if(family$family == "poisson" | family$family == "binomial"){

        copulaObject <- copula::normalCopula(copula::P2p(clusterCorrelationMatrix), dim = clusterSize, dispstr = "un")
                      
        eta <- mn.Y
                                         
        mu <- family$linkinv(eta)
        quantiles <- copula::rCopula(numClusters,copulaObject)
        quantiles <- as.vector(t(quantiles))

        if(family$family == "poisson" ){
            Y <- qpois(quantiles, mu)
        } else{
            Y <- qbinom(quantiles, size = 1, prob = mu)
        }

        
    } else {
        stop("Given family not yet implemented, please use either gaussian or poisson")
    }

    
    
      
    
    result <- list(y = Y,
                   x = X,
                   clusterID = rep(1:numClusters, each = clusterSize),
                   groupID = rep(1:numGroups, each = groupSize))

    if(!is.null(numMainEffects)){
        groupID = rep(1:((p^2 +p) /2), each = 1)
        result$xMainEff <- XmainEff
        result$interactionID <- interactionID
    }
    
    result
    
}
