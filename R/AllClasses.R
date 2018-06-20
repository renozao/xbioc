# Common S4 Classes
# 
# Author: Renaud Gaujoux
# Created: Mar 12, 2013
###############################################################################

#' @import methods
#' @import Biobase
#' @import AnnotationDbi
#' @import pkgmaker
#' @include utils.R
NULL

#' Generic Bioconductor Utilities
#' 
#' Utilities to test Bioconductor object classes.
#'
#' @name bioc-utils
NULL

#' @describeIn bioc-utils tells if an object inherits from the class [Biobase::ExpressionSet]. 
#' @param x object
#' @export
isExpressionSet <- function(x){
    is(x, 'ExpressionSet')
}

# union class for matrix-like data
setClassUnion('MatrixData', c('matrix', 'ExpressionSet'))

#' @describeIn bioc-utils tells if an object is a matrix-like data.
#' @export
isMatrixData <- function(x){
    is(x, 'MatrixData') || isExpressionSet(x)
}
