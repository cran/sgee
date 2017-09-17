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
#################################################################################' deltaFinder
#' 
#' Function to find the appropriate group update. A potential
#' update value is given and then the corresponding sparse group
#' lasso penalty value is calculated.
#' 
#' Internal function used by bisee.
#' 
#' @param a Possbile update value.
#' @param lambda1 Tuning parameter pertaining to importance of groups.
#' @param lambda2 Tuning parameter pertaining to importance of individuals.
#' @param gamma Scaling / thresholding value.
#' @return L2 Norm minus one of a given delta group defined by
#' given lambda1, lambda2, and gamma.
#' 
#' @note Internal function.
#' @author Gregory Vaughan
#' @keywords internal
deltaFinder <-
    function(a,
             lambda1,
             lambda2,
             gamma){
        return(sqrt(sum((deltaValue(a, lambda1, lambda2, gamma))^2))-1)
    }
