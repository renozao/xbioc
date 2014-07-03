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

# union class for matrix-like data
setClassUnion('MatrixData', c('matrix', 'ExpressionSet'))

# tests if an object is a matrix-like data.
isMatrixData <- function(x) is(x, 'MatrixData')
