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
    function(object, x){
      if( isString(x) ){
        if( !x %in% varLabels(object) )
          stop(sprintf("Could not find phenotypic variable '%s' in ExpressionSet object.", x))
        x <- object[[x]]
      }
      callGeneric(object = exprs(object), x = x)
      
    })




#' Melt an ExpressionSet for Use with ggplot 
#' 
#' Melt the expression data of an \code{\link[Biobase]{ExpressionSet}}, 
#' also binding any corresponding phenotypic data.
#' 
#' @param data an \code{ExpressionSet} object.
#' @param ... arguments passed to [reshape2::melt.array].
#' Named arguments, will result in separate variables in the returned \code{data.frame}.
#' @param pData phenotypic data to attach to the result \code{data.frame}.
#' @inheritParams reshape2::melt
#' @param extra list of other matrix-like data to melt and `rbind.fill` to the melted primary data.
#' 
#' @import reshape2
#' @importFrom plyr rbind.fill 
#' @export
melt.ExpressionSet <- function(data, ..., pData = pData(data), na.rm = FALSE, value.name = "value", extra = list()){
  
  # melt expression data
  pd <- NULL
  if( missing(pData) && hasMethod('pData', class(data)) ){
    pd <- getFunction('pData')(data)
  }else if( !missing(pData) ) pd <- pData
  
  X <- exprs(data)
  df <- melt(X, na.rm = na.rm, value.name = value.name, ...)
  
  if( length(extra) ){
    if( is.null(names(extra)) ) names(extra) <- rep('', length(extra))
    mextra <- ldply(seq_along(extra), function(i){
          n <- names(extra)[i]
          dx <- extra[[i]]
          if( ncol(dx) != ncol(data) ) 
            stop(sprintf("Invalid %i-th extra data [%s]: number of column [%i] does not match those from data [%i]", i, n, ncol(dx), ncol(data)))
          if( !nzchar(n) ) n <- colnames(df)[1]
          # force column names from data
          colnames(dx) <- colnames(data)
          # melt
          res <- melt(dx, na.rm = na.rm, value.name = value.name)
          # force same name for data column variable
          colnames(res)[1:2] <- c(n, colnames(df)[2])
          res
        })
    # add to result
    df <- rbind.fill(df, mextra)
  }
  
  # match pheno data
  if( !is.null(pd) ){
    i <- match(df[[2L]], colnames(data))
    stopifnot( !anyNA(i) )
    pd <- pd[i, , drop = FALSE]
    rownames(pd) <- NULL
    df <- cbind(df, pd)
  }
  
  df
  
}


#' Melting Assay Data
#' 
#' @param object an object with a suitable [Biobase::assayDataElement] method,
#' typically an [Biobase::ExpressionSet-class] object
#' @param elements an optional character vector that indicates the names of the 
#' assay data elements to include in the result
#' @param varnames a character vector with 2 unique elements that indicates the name of the columns
#' that will hold the feature and sample names respectively.
#' 
#' @return a `data.frame` object with the following structure:
#' 
#'   * character columns `feature_id` and `sample_id` containing the feature and sample names respectively;
#'   * other columns are numeric and correspond to each assay data element in the input object.
#' 
#' @importFrom checkmate assert_subset assert_character
#' @importFrom reshape2 dcast
#' @export
assayDataMelt <- function(object, elements = NULL, varnames = c("feature_id", "sample_id")){
  all_elements <- assayDataElementNames(object)
  if( is.null(elements) ) elements <- all_elements
  assert_subset(elements, all_elements)
  assert_character(varnames, len = 2L, min.chars = 1L, unique = TRUE)
  
  df <- ldply(setNames(nm = elements), .id = "element", function(elt){
    melt(assayDataElement(object, elt = elt), varnames = varnames, value.name = "value", as.is = TRUE)
    
  })
  reshape2::dcast(df, feature_id + sample_id ~ element, value.var = "value")
  
}
