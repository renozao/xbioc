# Project: xbioc
# 
# Author: Renaud Gaujoux
# Created: Sep 23, 2015
###############################################################################


#' Extracting Control Probes
#' 
#' @param object expression data as a \code{matrix} or \code{ExpressionSet} object.
#' 
#' @export
featureIsControl <- function(object){
    
    ids <- featureNames(object)
    
    # Affymetrix
    res <- grepl("^AFFX", ids)
    
    res
    
}
