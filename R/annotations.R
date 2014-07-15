# Project: xbioc
# 
# Functions for biological feature annotations
#
# Author: Renaud Gaujoux
# Created: Jul 3, 2014
###############################################################################

#' Handling Object Annotations
#' 
#' \code{hasAnnotation} tells if an object has some -- non empty -- attached annotation.
#' 
#' @param object an object
#' @param ... extra parameters (currently not used) 
#' 
#' @rdname annotation-utils
#' @export
hasAnnotation <- function(object, ...){
    !is.null( getAnnotation(object, ...) )
}
#' \code{getAnnotation} try extracting embedded annotations.
#' By default, it returns \code{NULL} if the object contains no annotation.
#'
#' @rdname annotation-utils
#' @export
getAnnotation <- function(object, ...){
    UseMethod('getAnnotation')
}

#' @param null logical that indicates if an empty character string should 
#' be return as \code{NULL}.
#' @rdname annotation-utils
#' @S3method getAnnotation default 
getAnnotation.default <- function(object, ..., null=TRUE){
    ann <- if( hasMethod('annotation', class(object)) ) annotation(object)
            else attr(object, 'annotation')
    if( null && isString(ann) && nchar(ann) == 0L ) ann <- NULL
    else if( !null && is.null(ann) ) ann <- '' 
    ann
}

#' \code{setAnnotation} embbeds an annotation string into an object 
#' (typically the name of an annotation package, e.g., \code{'hgu133plus2.db'}).
#' 
#' \code{setAnnotation} uses a suitable \code{annotation<-} method if it 
#' exists, or sets the annotation string into attribute \code{'annotation'}.
#' It returns the modified object.
#'
#' @param value new annotation string, e.g., \code{'hgu133plu2.db'}.
#'  
#' @rdname annotation-utils
#' @export
setAnnotation <- function(object, value, ...){
    UseMethod('setAnnotation')
}

#' @param as.attribute logical that indicates that the annotation string can be
#' stored as an attribute, if no suitable \code{`annotation <-`} method is found.
#'  
#' @rdname annotation-utils
#' @S3method setAnnotation default
setAnnotation.default <- function(object, value, as.attribute = TRUE, ...){
    
    if( length(extras <- list(...)) ) 
        stop("Unused arguments ", str_out(sapply(extras, class), Inf, use.names = TRUE))
    
    # pre-process
    if( is.null(value) ) value <- ''
    else if( length(value) > 1L ) value <- .glue_anndb(value)
    
    # set annotation
    if( hasMethod('annotation<-', c(class(object), class(value))) ) annotation(object) <- value
    else if( as.attribute ) attr(object, 'annotation') <- value
    else stop("Could not set annotation string: no suitable method `annotation <-` found for signature (", class(object), ",", class(value), ")")
    
    # return modified object
    object
}


#' @S3method getAnnotation list
getAnnotation.list <- function(object, ...){
    
    # handle list of ExpressionSet objects
    if( all(sapply(object, isExpressionSet)) ){
        gpl <- sapply(object, annotation)
        return( gpl )
    }
    NextMethod()
}


#' @S3method setAnnotation list
setAnnotation.list <- function(object, value, force = FALSE, ...){
    
    # handle list of ExpressionSet objects
    if( all(sapply(object, isExpressionSet)) ){
        gpl <- sapply(object, annotation)
        # try set annotation if provided
        if( all(gpl == gpl[1L]) || force ){
            sapply(seq_along(object), function(i){
                        annotation(object[[i]]) <<- value
                    })
        }else warning("Could not set provided annotation '", value,"': datasets are not from the same platform [", str_out(gpl, Inf), "]")
        return( object )
    }
    NextMethod()
}

#' Simple Feature Annotation 
#' 
#' @param x ENTREZ gene ids
#' @param annotation Annotation package to use.
#' It can be the name of an annotation package or of an organism supported by 
#' \code{\link{biocann_orgdb}}. 
#' @param link Indicates the type of online resource to link gene identifiers.
#' If not \code{'none'}, then ENTREZID and Symbols columns contain HTML links 
#' (<a> tags) that links to the corresponding gene page.
#' 
#' @seealso \code{\link{biocann_orgdb}}
#' 
#' @export
#' @examples 
#' 
#' # Entrez IDs
#' geneInfo(1:5)
#' geneInfo(1:5, link = 'bioGPS')
#' 
#' # probe ids
#' geneInfo(c("1007_s_at", "1053_at", "117_at", "121_at", "not_a_probe_id"), 'hgu133plus2.db')
#' geneInfo(c("1007_s_at", "1053_at", "117_at", "121_at", "not_a_probe_id"), 'hgu133plus2.db', link = 'bioGPS')
#' 
#' 
geneInfo <- function(x, annotation = 'human', link = c('none', 'bioGPS', 'NCBI')){
    
    use_org <- FALSE
    if( !is.annpkg(annotation) ){
        use_org <- TRUE
        db <- biocann_orgdb(annotation)
        annotation <- db$org.db
    }
    irequire(annotation, ptype = 'BioCann')
    
    x <- as.character(x)
    ez <- if( use_org ) x else bimap_lookup(x, biocann_object('ENTREZID', annotation), multiple = FALSE)
    symb <- bimap_lookup(x, biocann_object('SYMBOL', annotation), multiple = FALSE)
    ez[is.na(symb)] <- NA
    desc <- bimap_lookup(x, biocann_object('GENENAME', annotation), multiple = FALSE)
    
    df <- data.frame(ENTREZID = ez, Symbol = symb, Description = desc, stringsAsFactors = FALSE)
    rownames(df) <- x
    
    # convert to HTML links
    link <- match.arg(link)
    if( link != 'none' ){
        tolink <- c('ENTREZID', 'Symbol')
        base <- switch(link
                , bioGPS = "http://biogps.org/gene/"
                , NCBI = "http://www.ncbi.nlm.nih.gov/gene/")
        .link <- function(x, id = x){
            ok <- !is.na(id) & nzchar(x)
            x[ok] <- sprintf('<a target="_" href="%s%s">%s</a>', link, base, id[ok], x[ok])
            x
        }
        df[tolink] <- sapply(df[tolink], .link, id = df$ENTREZID, simplify = FALSE)
    }
    
    df
}
