# Common set of S4 generics
# 
# Author: Renaud Gaujoux
# Created: 28 Nov 2012
###############################################################################

#' @include AllClasses.R
NULL

#' Extracting Feature Names
#' 
#' The \pkg{xbioc} package provides extra methods for the generic 
#' \code{\link[Biobase]{featureNames}} and \code{\link[Biobase]{sampleNames}}, 
#' that extend the original Bioconductor interface.
#' 
#' @inheritParams Biobase::featureNames
#' 
#' @rdname Bioc-generics
#' @export
#' @inline
setGeneric('featureNames', package='Biobase')

#' @rdname Bioc-generics
#' @export
#' @inline
setGeneric('featureNames<-', package='Biobase')

#' Returns the row names of \code{object}.
setMethod('featureNames', 'matrix', function(object) rownames(object) )

#' Set the row names of \code{object}.
setReplaceMethod('featureNames', 'matrix', 
   function(object, value){
        rownames(object) <- value
        object
   }
)

#' @rdname Bioc-generics
#' @export
#' @inline
setGeneric('sampleNames', package='Biobase')

#' @rdname Bioc-generics
#' @export
#' @inline
setGeneric('sampleNames<-', package='Biobase')

#' Set the row names of \code{object}.
#' @aliases sampleNames,matrix,ANY
setReplaceMethod('sampleNames', 'matrix', 
    function(object, value){
        colnames(object) <- value
        object
    }
)


#' Returns the column names of \code{object}.
setMethod('sampleNames', 'matrix', function(object) colnames(object) )

#' @export
#' @rdname Bioc-generics
setGeneric('exprs', package = 'Biobase')

#' Simply returns \code{object}.
#' This method is defined so the generic can be called on both matrix and 
#' \code{ExpressionSet} objects.
#' 
#' @rdname Bioc-generics
setMethod('exprs', 'matrix', function(object) object )

#' @export
#' @rdname Bioc-generics
setGeneric('exprs<-', package = 'Biobase')

#' Simply assigns `value` to `object`.
#' This method is defined so the generic can be called on both matrix and 
#' \code{ExpressionSet} objects.
#' 
#' @rdname Bioc-generics
setReplaceMethod('exprs', 'matrix', function(object, value){
    object <- value
    object
  })


#' Enhanced Subsetting for Matrix-like Data 
#' 
#' These methods subset Matrix-like data only keeping rows 
#' whose names are in common with some other data.
#' 
#' @inheritParams base::intersect
#' @name intersect
NULL

#' This method is equivalent to \code{\link[base]{subset}(x, y)}.
#' @rdname intersect
#' @export 
setMethod('intersect', signature(x='MatrixData', y='logical'), 
        function(x, y){
            subset(x, y)
        }
)
#' Subset a matrix-like object by only keeping the rows whose 
#' names are in a given reference character vector.
#' @rdname intersect
#' @export 
setMethod('intersect', signature(x='MatrixData', y='character'), 
        function(x, y){
            intersect(x, featureNames(x) %in% y)
        }
)
#' Subset a matrix-like object by only keeping the rows whose 
#' names are in common with another matrix-like data.
#' 
#' This is a shortcut for \code{intersect(x, featureNames(y), ...)}.
#' @rdname intersect
#' @export 
setMethod('intersect', signature(x='MatrixData', y='MatrixData'), 
        function(x, y){
            intersect(x, featureNames(y))
        }
)
