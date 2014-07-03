# Project: xbioc
#
# Functions for gene expression data
#
# Author: Renaud Gaujoux
# Created: Jul 3, 2014
###############################################################################


#' Compare Features Sets
#'
#' Compares the set of features associated with an object, e.g., \code{\link[Biobase]{ExpressionSet}} 
#' or \code{\link[AnnotationDbi]{ChipDb}} objects.
#' 
#' @param x, y objects with associated featureNames
#' @param ... extra arguments eventually passed to \code{compareFeatures,character,character}. 
#' 
#' @export
#' @inline
setGeneric('compareFeatures', function(x, y, ...){ standardGeneric('compareFeatures') })

#' Workhorse function comparing the feature sets \code{x} and \code{y}. 
setMethod('compareFeatures', signature(x = 'character', y = 'character'), 
    function(x, y){
        common <- intersect(x, y)
        x.only <- setdiff(x, y)
        y.only <- setdiff(y, x)
        xl <- length(x)
        yl <- length(y)
        cat(sprintf("  x: %i\n  y: %i\n", xl, yl))
        cat(sprintf("  Common: %.2f%% x | %.2f%% y > %s\n", length(common)/ xl * 100, length(common)/ yl * 100
                , str_out(common, total = TRUE)))
        cat(sprintf("  x only: %.2f%% > %s\n", length(x.only)/ xl * 100, str_out(x.only, total = TRUE)))
        cat(sprintf("  y only: %.2f%% > %s\n", length(y.only)/ yl * 100, str_out(y.only, total = TRUE)))
})

setMethod('compareFeatures', signature(x = 'ANY', y = 'ANY'), 
    function(x, y, ...){
        if( !is.character(x) ) x <- featureNames(x)
        if( !is.character(y) ) y <- featureNames(y)
        callGeneric()
})

setMethod('compareFeatures', signature(x = 'matrix'), 
    function(x, y, ...){
        x <- rownames(x)
        callGeneric()
})

setMethod('compareFeatures', signature(y = 'matrix'), 
    function(x, y, ...){
        y <- rownames(y)
        callGeneric()
})

setMethod('compareFeatures', signature(x = 'ChipDb'), 
    function(x, y, ...){
        x <- keys(x)
        callGeneric()
})

setMethod('compareFeatures', signature(y = 'ChipDb'), 
    function(x, y, ...){
        y <- keys(y)
        callGeneric()
})

