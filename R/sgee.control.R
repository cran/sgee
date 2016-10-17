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
#'
#' @return A list containing all of the parameter values.
#' 
#' @author Gregory Vaughan
#' @export
sgee.control <- function(maxIt = 200 , epsilon = 0.05,
                         stoppingThreshold =  NULL){

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

    list(maxIt = maxIt,
         epsilon = epsilon,
         stoppingThreshold = stoppingThreshold)
}


