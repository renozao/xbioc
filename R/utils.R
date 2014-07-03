# Project: GEOdb
# 
# Author: renaud
# Created: Nov 18, 2013
###############################################################################

#' @include AllClasses.R
NULL

.glue_anndb <- function(x) paste0(x, collapse=', ')


#' Detect Log-transformed Data
#' 
#' \code{is_logscale} tells if some numeric data is in log scale, 
#' e.g., normalized microarray data, using the same heuristic as GEO2R.
#' 
#' The data needs to be of reasonable size and variance for the detection 
#' heuristic to work correctly.
#' 
#' @param x a numeric data object (matrix, vector, ExpressionSet) 
#' 
#' @source \url{www.ncbi.org/geo}
#' 
#' @export
#' @examples
#' 
#' x <- matrix(rnorm(20*10, mean=500), 20, 10)
#' is_logscale(x)
#' is_logscale(log_transform(x))
#' 
is_logscale <- function(x){
    
    ex <- if( isExpressionSet(x) ) exprs(x) else x
    # check log2 transform
    #ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    !LogC
#	if (LogC) { ex[which(ex <= 0)] <- NaN
#		exprs(gset) <- log2(ex) }
}

#' \code{log_transform} apply a log transformation to the data.
#' Negative values get assigned the value \code{\link{NaN}}.
#' 
#' @param base log base to use.
#' 
#' @export
#' @rdname is_logscale
log_transform <- function(x, base=2){
    
    ex <- if( isExpressionSet(x) ) exprs(x) else x
    
    # transform
    ex[which(ex <= 0)] <- NaN
    ex <- log(ex, base)
    
    # return same type of object
    if( isExpressionSet(x) ) exprs(x) <- ex
    else x <- ex
    x
}

#' Generic Bioconductor Utilities
#' 
#' \code{isExpressionSet} tells if an object inherits from the class \code{\link[Biobase]{ExpressionSet}}.
#' 
#' @param x object 
#' @export
isExpressionSet <- function(x){
    is(x, 'ExpressionSet')
}

