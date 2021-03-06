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
#################################################################################' sgee: Stagewise Generalized Estimating Equations 
#' 
#' Provides functions to perform Boosting / Functional Gradient Descent /
#' Forward Stagewise regression with grouped covariates setting using
#' Generalized Estimating Equations.
#' 
#' \tabular{ll}{ Package: \tab sgee\cr Type: \tab Package\cr
#' Version: \tab 0.6-0\cr Date: \tab 2018-01-08\cr License: \tab GPL (>= 3)
#' \cr } sgee provides several stagewise regression approaches
#' that are designed to address variable selection with grouped covariates
#' in the context of
#' Generalized Estimating Equations. Given a response and design matrix
#' stagewise techniques perform a sequence of small learning steps
#' wherein a subset of the covariates are selected as being the
#' most important at that iteration and are then subsequently updated
#' by a small amount, epsilon. different techniques this optimal update
#' in different ways that achieve different structural goals (i.e.
#' groups of covariates are fully included or not).
#'
#' The resulting path can then be analyzed to determine an optimal
#' model along the path of coefficient estimates. The
#' \code{analyzeCoefficientPath} function provides such
#' functionality based on various
#' possible metrics, primarily focused on the Mean Squared Error.
#' Furthermore, the \code{plot.sgee} function can be used to examine the
#' path of coefficient estimates versus the iteration number, or some
#' desired penalty.
#' 
#' @name sgee-package
#' @aliases sgee-package sgee
#' @docType package
#' @author Gregory Vaughan [aut, cre],
#' Kun Chen [ctb], Jun Yan [ctb]
#' 
#' Maintainer: Gregory Vaughan <gregory.vaughan@uconn.edu>
#' @references Vaughan, G., Aseltine, R., Chen, K., Yan, J., (2017). Stagewise
#' Generalized Estimating Equations with Grouped Variables. Biometrics 73,
#' 1332-1342. URL: http://dx.doi.org/10.1111/biom.12669,
#' doi:10.1111/biom.12669.
#'
#' Vaughan, G., Aseltine, R., Chen, K., Yan, J., (2017). Efficient
#' interaction selection for clustered data via stagewise generalized
#' estimating equations.  Department of Statistics, University of
#' Connecticut. Technical Report.
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
#' 
#' Hastie, T., Tibshirani, R., and Friedman, J. (2009). The Elements
#' of Statistical Learning: Data Mining, Inference, and Prediction.
#' Springer, New York.
#' 
#' Liang, K.-Y. and Zeger, S. L. (1986). Longitudinal data analysis
#' using generalized linear models. Biometrika 73, 13--22.
#' @examples
#' 
#' 
#' 
#' #####################
#' ## Generate test data
#' #####################
#' 
#' ## Initialize covariate values
#' p <- 50 
#' beta <- c(rep(2.4,5),
#'           c(1.2, 0, 1.6, 0, .4),
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
#'                          intercept = 0)
#' 
#' 
#' coefMat1 <- hisee(y = generatedData$y, x = generatedData$x,
#'                   family = gaussian(),
#'                   clusterID = generatedData$clusterID,
#'                   groupID = generatedData$groupID, 
#'                   corstr="exchangeable",
#'                   control = sgee.control(maxIt = 100, epsilon = 0.2))
#' 
#' ## interceptLimit allows for compatibility with older R versions
#' coefMat2 <- bisee(y = generatedData$y, x = generatedData$x,
#'                   family = gaussian(),
#'                   clusterID = generatedData$clusterID,
#'                   groupID = generatedData$groupID, 
#'                   corstr="exchangeable", 
#'                   control = sgee.control(maxIt = 100, epsilon = 0.2,
#'                                          interceptLimit = 10),
#'                   lambda1 = .5,
#'                   lambda2 = .5)
#' 
#' 
#' par(mfrow = c(2,1))
#' plot(coefMat1)
#' plot(coefMat2)
#' 
NULL



