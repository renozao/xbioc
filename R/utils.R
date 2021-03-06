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
#' @param robust logical that indicates if the decision should be 
#' done robustly by removing the more extreme values.
#' @param extremes a length-2 numeric vector that gives the lower and upper 
#' quantiles used to remove extreme values. 
#' Only relevant when `robust = TRUE`.
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
is_logscale <- function(x, robust = TRUE, extremes = c(0.02, 0.98)){
    
    ex <- if( isExpressionSet(x) ) exprs(x) else x
    # check log2 transform
    #ex <- exprs(gset)

    # move negative values to positive
    if (any(ex < 0, na.rm = TRUE)) ex <- ex - min(ex, na.rm = TRUE)
    
    # remove extreme values 
    if( robust ){
      stopifnot( length(extremes) == 2L )
      q_extreme <- quantile(ex, probs = sort(extremes), na.rm = TRUE)
      ex[ex < q_extreme[1L]] <- NA_real_
      ex[ex > q_extreme[2L]] <- NA_real_
      
    }
    ##

    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
    LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    !LogC
#	if (LogC) { ex[which(ex <= 0)] <- NaN
#		exprs(gset) <- log2(ex) }
}

#' @describeIn is_logscale tests if there are outliers to the log-scale distribution.
#' @param ... other arguments passed to [is_logscale], like argument `extremes`
#' to control the extreme values that are removed.
#' @export 
has_logscale_outliers <- function(x, ...){
  is_logscale(x, robust = TRUE, ...) != is_logscale(x, robust = FALSE, ...)
  
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

#' Initialising a Mapping List
#' 
#' @param x object used as right keys
#' @return a list with names \code{x} mapping to NA values.
#' 
#' @export
#' 
#' @examples 
#' NAmap(letters[1:5])
#' 
NAmap <- function(x){
    as.list(setNames(rep(NA, length(x)), x))
}


#' Appending Columns to Data Frames
#' 
#' This function performs on \code{data.frame} objects a similar operation as 
#' \code{\link[base]{append}} does on \code{list} objects.
#' That is it enables prepending/inserting/appending columns at specific positions.
#' 
#' @param x a \code{data.frame} object
#' @param ... variables to be appended to \var{x}.
#' @param after index or column name after which the new variables should be inserted.
#' Default is to append the columns at the end of the \code{data.frame}.
#' @param names optional names for the newly inserted columns.
#' @param stringsAsFactors logical that indicates if the newly inserted string variables should be converted to 
#' factors.
#' 
#' @return a `data.frame`
#' 
#' @export
#' @examples 
#' 
#' x <- data.frame(a = 1:4, b = letters[1:4])
#' # append at the end
#' df_append(x, c = runif(4))
#' 
#' # append at the beginning
#' df_append(x, c = runif(4), after = 0L)
#' 
#' # append after first column
#' df_append(x, c = runif(4), after = 1L)
#' df_append(x, c = runif(4), after = 'a')
#' 
df_append <- function(x, ..., after = length(x), names = NULL, stringsAsFactors = getOption('stringsAsFactors', FALSE)){
  
  addon <- data.frame(..., stringsAsFactors = stringsAsFactors)
  if( isString(after) ) after <- which(colnames(x) == after)
  if( !length(after) ){
    warning("Argument 'after' has zero-length: value will be concatenated after last column.")
    after <- length(x)
  }
  n <- length(x)
  res <- if( after <= 0 ) cbind(addon, x)
      else if( after >= n ) cbind(x, addon)
      else cbind(x[seq(after)], addon, x[-seq(after)])
  
  # apply colnames if requested
  if( !is.null(names) ){
    colnames(res)[seq(after+1, after+length(addon))] <- names
  }
  
  res
  
}

#' Computing a Short Digest String
#' 
#' @param x an object
#' @param n a single integer/numeric that indicates the desired length of the digest string.
#' Use `n = Inf` to obtain full digest strings.
#' @param ... other arguments passed to [digest::digest]
#' 
#' @return a character vector of digest strings.
#' 
#' @importFrom digest digest
#' @export 
#' @examples
#' 
#' x <- list(1, 2, 3)
#' digestN(x)
#' # longer diest
#' digestN(x, 15)
#' # full digest
#' digestN(x, Inf)
#' 
digestN <- function(x, n = 7, ...){
  d <- digest(x, ...)
  if( is.finite(n) ) d <- substr(d, 1, n)
  d
  
}
