# Project: xbioc
#
# Flow Cytometry Utils 
#
# Author: Renaud Gaujoux
# Created: Oct 13, 2015
###############################################################################


#' Change Reference Populations in Flowcytometry Data
#' 
#' @param x flow data
#' @param ref indicates the reference population, that is either \code{'parent'} or
#' \code{'total'}
#' @export
`flow_reference<-` <- function(x, ref, parent, ...){
    UseMethod('flow_reference<-')
}

#' @export 
`flow_reference<-.default` <- function(x, parent, ..., value){
   
    # match ref mode
    ref <- value
    if( !is.null(names(parent)) ){
        rownames(x) <- names(parent)
    }
    
    subsets <- rownames(x)
    stopifnot( !anyDuplicated(subsets) )
    
    if( ref %in% subsets ){
        top_i <- which(subsets == ref)
        ref <- 'total'
           
    }else top_i <- which(is.na(parent))
    
    ref <- match.arg(ref, c('parent', 'total'))
    pct_scale <- 1
    if( sum(x>1) / length(x) > .5 ) pct_scale <- 100
    message("* Scaling factor: ", pct_scale)
    message("* Computing proportions to ", ref)
    
    res <- switch(ref,
        total = { # compute proportions to reference population
            flowTot <- x / pct_scale
            # force missing value for subsets with missing parent
            flowTot[setdiff(which(!parent %in% subsets), top_i), ] <- NA
            
            # start from top population
            ct <- subsets[top_i]
            while( length(i <- which(parent %in% ct)) ){ # look up cell types with parent in current level 
                  ip <- match(parent[i], subsets)
                  stopifnot( !anyNA(ip) )
                  # compute proportion relative to parent
                  flowTot[i, ] <- flowTot[i, ] * flowTot[ip, ]
                  # process to next level: children of current subsets
                  ct <- subsets[i]
                  
            }
            flowTot
        }
    , parent = {
        x
    })

    # return result
    res
    
}

#' @export 
`flow_reference<-.ExpressionSet` <- function(x, parent = x[['parentID']], ..., value){
   ref <- value
   if( isString(ref) ) ref <- x[[ref]] 
   x <- exprs(x)
   NextMethod()
}

#' Read All FCS Files in a Directory
#' 
#' Read all the FCS files in a directory, using 
#' \code{\link[flowCore]{read.FCS}}.
#' 
#' @param path path to the directory to load
#' @param ... arguments passed to \code{\link[flowCore]{read.FCS}}
#' @param pattern file pattern to match.
#' Files are matched using regular expressions via 
#' \code{\link{list.files}} in case insensitive mode.
#' @param recursive logical that indicates if sub-directory should 
#' be looked up.
#' 
#' @export
read.fcs <- function(path, ..., pattern = "\\.fcs$", recursive = FALSE){
    qlibrary('flowCore')
    sapply(list.files(path, pattern = pattern, full.names = TRUE
                        , ignore.case = TRUE, recursive = recursive)
            , flowCore::read.FCS, ..., simplify = FALSE)
    
}

