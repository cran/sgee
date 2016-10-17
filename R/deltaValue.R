
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
#################################################################################' deltaValue
#' 
#' Function to calculate an update under a given configuration.
#' 
#' Internal function used by bisee.
#' 
#' @param a Subvector of Estimating Equations that pertains to a particular
#' group.
#' @param lambda1 Tuning parameter pertaining to importance of groups.
#' @param lambda2 Tuning parameter pertaining to importance of individuals.
#' @param gamma Scaling / thresholding value.
#' @return The form of a group update for given lambda1, lambda2, and gamma. 
#' 
#' @note Internal function.
#' @author Gregory Vaughan
#' @keywords internal 
deltaValue <-
    function(a,
             lambda1,
             lambda2,
             gamma){
        
        
        result <- sign(a)*pmax(abs(a)-gamma*lambda2, 0)/(gamma*lambda1*sqrt(length(a)))
        return(result)
    }
