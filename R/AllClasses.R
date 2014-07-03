# Common S4 Classes
# 
# Author: Renaud Gaujoux
# Created: Mar 12, 2013
###############################################################################

#' @import methods
#' @import Biobase
#' @import AnnotationDbi
#' @import pkgmaker
NULL

#' Generic Bioconductor Utilities
#' 
#' \code{isExpressionSet} tells if an object inherits from the class \code{\link[Biobase]{ExpressionSet}}.
#' 
#' @param x object 
#' @export
isExpressionSet <- function(x){
    is(x, 'ExpressionSet')
}

# union class for matrix-like data
setClassUnion('MatrixData', c('matrix', 'ExpressionSet'))

#' \code{isExpressionSet} tells if an object is a matrix-like data.
#' @export
#' @rdname isExpressionSet
isMatrixData <- function(x){
    is(x, 'MatrixData') || isExpressionSet(x)
}
