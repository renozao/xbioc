# Utility for ExpressionSet objects
# 
# Author: Renaud Gaujoux
###############################################################################


#' Retrieve Phenotypic Variable from Data
#' 
#' This is a generic function to extract a phenotypic variable from data objects.
#' It aims is to encompass variables that are directly defined as vectors or,  
#' for example, specified as a name that is looked up in the `phenoData` variables
#' of an `ExpressionSet` object.
#' 
#' @param object an `ExpressionSet` object.
#' @param x variable as an atomic vector or as a string corresponding to the 
#' name of a phenotypic variable to be looked up in `phenData(object)`.
#' 
#' @export
setGeneric('pVar', function(object, x, ...) standardGeneric('pVar') )

#' @describeIn pVar simply returns `x`, after checking that it has the same length 
#' as the number of columns in `x`.
setMethod('pVar', signature(object = 'matrix'), 
    function(object, x){    
      if( length(x) != ncol(object) )
        stop(sprintf("Invalid phenotypic variable length [%s]: must be of length the number of samples in data [%s].", length(x), ncol(object)))
      x
    })
#' @describeIn pVar either return `x` or the phenotypic variable `pData(object)[[x]]`
#' if it exists.
setMethod('pVar', signature(object = 'ExpressionSet'), 
    function(object, x, check = TRUE, ...){
      if( isString(x) ){
        if( !x %in% varLabels(object) )
          stop(sprintf("Could not find phenotypic variable '%s' in ExpressionSet object.", x))
        x <- object[[x]]
      }
      callGeneric(object = exprs(object), x = x, ...)
      
    })



