% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.sgeeSummary.R
\name{print.sgeeSummary}
\alias{print.sgeeSummary}
\title{\code{print} function for sgee summaries}
\usage{
\method{print}{sgeeSummary}(x, ...)
}
\arguments{
\item{x}{An object of the \code{sgeeSummary} class, produced by
applying the \code{summary} function to an object of class \code{sgee}.}

\item{...}{Not currently used}
}
\description{
Provides implementation of \code{print} function for summaries of sgee
objects.
}
\examples{

#####################
## Generate test data
#####################

## Initialize covariate values
p <- 50 
beta <- c(rep(2,5),
          c(1, 0, 1.5, 0, .5),
          rep(0.5,5),
          rep(0,p-15))
groupSize <- 5
numGroups <- length(beta)/groupSize


generatedData <- genData(numClusters = 50,
                         clusterSize = 4,
                         clusterRho = 0.6,
                         clusterCorstr = "exchangeable",
                         yVariance = 1,
                         xVariance = 1,
                         numGroups = numGroups,
                         groupSize = groupSize,
                         groupRho = 0.3,
                         beta = beta,
                         family = gaussian(),
                         intercept = 1)

genDF <- data.frame(generatedData$y, generatedData$x)
names(genDF) <- c("Y", paste0("Cov", 1:p))
coefMat <- hisee(formula(genDF), data = genDF,
                 family = gaussian(),
                 clusterID = generatedData$clusterID,
                 groupID = generatedData$groupID, 
                 corstr="exchangeable", 
                 maxIt = 50,
                 epsilon = .5)

sgeeSum <- summary(coefMat)
print(sgeeSum)

}
\author{
Gregory Vaughan
}
