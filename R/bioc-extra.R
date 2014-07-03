# Bioconductor extensions/utility functions
# 
# Author: Renaud Gaujoux
# Creation: 10 Jan 2012
###############################################################################

#' @include AllGenerics.R 
NULL

#' Numeric Computations on ExpressionSet objects 
#' 
#' @description
#' The \pkg{CellMix} defines some generics and methods to apply numeric transformations 
#' to \code{ExpressionSet} objects, which is convenient when working on gene expression 
#' deconvolution algorithms, where scale (log/linear) may matter. 
#'  
#' \code{log} log-transforms the expression matrix of \code{\link{ExpressionSet}} objects.
#' 
#' @param x an \code{ExpressionSet} object.
#' @param ... extra arguments passed to subsequent calls, usually of the corresponding  
#' method in the \pkg{stats} package.
#' 
#' @rdname ExpressionSet-compute
#' @export
setMethod('log', 'ExpressionSet'
        , function(x, ...){
            exprs(x) <- log_transform(exprs(x), ...)
            x
        }
)

#' \code{expb} applies an entry-wise exponential transformation to the expression matrix 
#' of \code{\link{ExpressionSet}}, in a specific base.
#' 
#' @export
#' @rdname ExpressionSet-compute
#' @inline
setGeneric('expb', function(x, ...) standardGeneric('expb') )
#' @param base log base to use.
#' @export
setMethod('expb', 'numeric'
        , function(x, base=exp(1)){	
            if( !missing(base) ) exp(x * log(base)) else exp(x)
        }
)
#' @export
setMethod('expb', 'matrix', selectMethod('expb', 'numeric'))
#' @export
setMethod('expb', 'ExpressionSet'
        , function(x, ...){	
            exprs(x) <- expb(exprs(x), ...) 
            x
        }
)

#' \code{exp} is equivalent to \code{expb(x, exp(1))}. 
#' 
#' @export
#' @rdname ExpressionSet-compute
setMethod('exp', 'ExpressionSet', function(x) expb(x) )

#' \code{range} computes the range of expression values from 
#' \code{\link{ExpressionSet}} objects.
#' 
#' @export
#' @rdname ExpressionSet-compute
setMethod('range', 'ExpressionSet'
        , function(x, ..., na.rm = FALSE){
            range(exprs(x), ..., na.rm=na.rm)
        }
)

#' \code{quantile} computes the range of expression values in 
#' \code{\link{ExpressionSet}} objects.
#' 
#' @S3method quantile ExpressionSet
#' @rdname ExpressionSet-compute
quantile.ExpressionSet <- function(x, ...){
    quantile(exprs(x), ...)
}

#' Combining Expression Matrices
#' 
#' The method \code{cbind.ExpressionSet} combines sample expression 
#' profiles from multiple \code{\link{ExpressionSet}} or matrix objects.
#' 
#' The expression matrices must be exactly of the same dimensions.
#' For the end result to be meaningful, one probably wants the row names 
#' to match as well, i.e. both object contain the same features, in the same order. 
#' However no check is done for this.
#' 
#' Note that the returned \code{ExpressionSet} object has no sample or feature 
#' annotations.
#' 
#' @param ... series of \code{ExpressionSet} and/or matrix objects.
#' @inheritParams base::cbind
#' 
#' @return an \code{ExpressionSet} object
#' 
#' @S3method cbind ExpressionSet
#' @export
cbind.ExpressionSet <- function(..., deparse.level = 1){
    objects <- list(...)
    nobjects <- length(objects)
    if( nobjects == 1L && is.list(objects[[1L]]) ){
        objects <- objects[[1L]]
        nobjects <- length(objects)
    }
    if( !nobjects ) stop("Could not generate ExpressionSet: empty input list.")
    
    out <- objects[[1]]
    if( !isExpressionSet(out) )
        stop("Could not generate ExpressionSet: invalid first element [", class(out), "]")
    #other <- names(objects[[1]]$other)
    # TODO: check feature names
    # init pheno data
    pd <- pData(out)
    out <- exprs(out)
    pd_class <- sapply(pd, class)
    pd <- sapply(pd, as.character)
    if (nobjects > 1){ 
        lapply(2:nobjects, function(i){
                    # expression values
                    o <- objects[[i]]
                    out <<- cbind(out, if( is(o, 'ExpressionSet') ) exprs(o) else o )
                    # binds pheno data
                    pd_o <- sapply(pData(o), as.character)
                    pd <<- rbind(pd, pd_o)
                    NULL
                })
    }
    #pd <- mapply(as, pd, pd_class)
    sn <- unlist(lapply(objects, sampleNames), use.names = FALSE)
    pd <- as.data.frame(pd, row.names = sn)
#    str(pd)
    ExpressionSet(out, phenoData = AnnotatedDataFrame(pd), annotation = annotation(objects[[1]]))
}


#' @S3method droplevels ExpressionSet
droplevels.ExpressionSet <- function(x, ...){
    # drop levels from annotation data
    pData(x) <- droplevels(pData(x))
    fData(x) <- droplevels(fData(x))
    
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

