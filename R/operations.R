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
#' The \pkg{xbioc} package defines some generics and methods to apply numeric transformations 
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
#' @importFrom stats quantile
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
#' @param .id name of the phenotypic variable created to old each dataset ID.
#' If argument are named then these are used as values.  
#' @inheritParams base::cbind
#' 
#' @return an \code{ExpressionSet} object
#' 
#' @S3method cbind ExpressionSet
#' @importFrom plyr rbind.fill.matrix
#' @export
cbind.ExpressionSet <- function(..., .id = '.id', deparse.level = 1){
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
    # check feature names
    fn <- sapply0(objects, featureNames)
    stopifnot(all(sapply(fn , identical, fn[[1L]])))
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
                    pd <<- rbind.fill.matrix(pd, pd_o)
                    NULL
                })
    }
    #pd <- mapply(as, pd, pd_class)
    sn <- unlist(lapply(objects, sampleNames), use.names = FALSE)
    pd <- as.data.frame(pd, row.names = sn)
    # add .id
    if( is.null(ids <- names(objects)) ){
        ids <- as.character(seq_along(objects))
    } 
    pd[[.id]] <- factor(unname(unlist(mapply(rep, ids, sapply(objects, ncol)))))
    # fData
    fd <- featureData(objects[[1L]])
#    rownames(fd) <- rownames()
#    str(pd)
    ExpressionSet(out, phenoData = AnnotatedDataFrame(pd), featureData = fd, annotation = annotation(objects[[1]]))
}


#' @S3method droplevels ExpressionSet
droplevels.ExpressionSet <- function(x, ...){
    # drop levels from annotation data
    pData(x) <- droplevels(pData(x))
    fData(x) <- droplevels(fData(x))
    
    x
}

#' @S3method rbind ExpressionSet 
#' @export
rbind.ExpressionSet <- function(..., .id = 'dataType'){
    
    DATA <- list(...)
    if( is.null(names(DATA)) ) names(DATA) <- paste0('dataset', seq_along(DATA))
    if( anyDuplicated(names(DATA)) )
        stop("Invalid expression data list: names should be unique.")
    qlibrary(plyr)
    # check features do not overlap
    ok <- lapply(names(DATA), function(n){
            lapply(setdiff(names(DATA), n), function(p){
                ok <- length(intersect(featureNames(DATA[[n]]), featureNames(DATA[[p]])))
                setNames(ok, paste0(n, '-', p))
            })
        })
    ok <- unlist(ok)
    if( any(ok>0) )
        stop("Overlapping expression data: datasets should be disjoint [", str_out(ok[ok>0], Inf, use.names = TRUE), "]")
    
    # temporary auxiliary column
    nameCol <- basename(tempfile('NAMES_'))
    
    # extract data
    mat <- as.matrix(ldply(DATA, function(x) exprs(x), .id = NULL))
    rownames(mat) <- unlist(lapply(DATA, featureNames))
    
    # extract and merge pheno data
    pd <- ldply(DATA, function(x){
        pd <- pData(x)
        pd[[nameCol]] <- sampleNames(x)
        pd
    }, .id = NULL)
    pd <- pd[!duplicated(pd[[nameCol]]), ]
    rownames(pd) <- pd[[nameCol]]
    pd <- pd[setdiff(names(pd), nameCol)]
    
    # extract and bind fetaure data
    fd <- ldply(DATA, function(x){
        fd <- fData(x)
        fd[[.id]] <- NULL
        fd[[nameCol]] <- featureNames(x)
        fd
    }, .id = .id)
    stopifnot( !anyDuplicated(fd[[nameCol]]) )
    rownames(fd) <- fd[[nameCol]]
    fd <- fd[setdiff(names(fd), nameCol)]
    
    # build object
    ExpressionSet(mat, phenoData = AnnotatedDataFrame(pd[colnames(mat), ])
            , featureData = AnnotatedDataFrame(fd[rownames(mat), ])
        )
    
}
