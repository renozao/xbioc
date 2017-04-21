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
#' @param x feature identifiers to annotate
#' @param annotation Annotation package to use.
#' It can be the name of an annotation package or of an organism supported by 
#' \code{\link{biocann_orgdb}}. 
#' 
#' It can also be a `AnnDbBimap` object that converts the features into `ENTREZID`, e.g., `org.Hs.egSYMBOL2EG`.
#' @param extras Indicates the type of information/resource to add for each feature.
#' Out-links to online resource can be added the prefix \code{'~'} to selected resources.
#' 
#' Note that the resources \code{'NCBI'} and \code{'bioGPS'} are \emph{de facto} online resources, 
#' and make ENTREZID and Symbol as links to the respective gene page. 
#' 
#' @seealso \code{\link{biocann_orgdb}}
#' 
#' @export
#' @examples 
#' 
#' # Entrez IDs
#' geneInfo(1:5)
#' geneInfo(1:5, extras = 'bioGPS')
#' 
#' # probe ids
#' ids <- c("1007_s_at", "1053_at", "117_at", "121_at", "not_a_probe_id")
#' geneInfo(ids, 'hgu133plus2.db')
#' geneInfo(ids, 'hgu133plus2.db', extras = 'bioGPS')
#' geneInfo(ids, 'hgu133plus2.db', extras = 'pathwayID')
#' 
#' 
geneInfo <- function(x, annotation = 'human', extras = c('biogps', 'ncbi', 'pathway', 'kegg')){
    
    # handle ExpressionSet objects
    if( isExpressionSet(x) ){
        if( missing(annotation) || is.null(annotation) ) annotation <- annotation(x)
        x <- fData(x)
    }
    
    # use rownames from data.frame and append/cbind original data at the end 
    bind_data <- NULL
    if( is.data.frame(x) ){
      bind_data <- x
      x <- rownames(x)
    }
    
    x0 <- x
    use_org <- FALSE
    add_id <- TRUE
    if( !is.annpkg(annotation) ){
        use_org <- TRUE
        # assume that it will convert to ENTREZID
        if( is(annotation, 'AnnDbBimap') ){
          x <- bimap_lookup(as.character(x), annotation, multiple = FALSE)
        } else add_id <- FALSE
        # use annotation from organism (key: ENTREZID)
        annotation <- biocann_orgdb_pkgname(annotation)
        
    }
    irequire(annotation, ptype = 'BioCann')
    
    x <- as.character(x)
    ez <- if( use_org ) x else bimap_lookup(x, biocann_object('ENTREZID', annotation), multiple = FALSE)
    symb <- bimap_lookup(x, biocann_object('SYMBOL', annotation), multiple = FALSE)
    ez[is.na(symb)] <- NA
    desc <- bimap_lookup(x, biocann_object('GENENAME', annotation), multiple = FALSE)
    df <- data.frame(ENTREZID = ez, SYMBOL = symb, Description = desc, stringsAsFactors = FALSE)
    if( add_id ) df$ID <- x0
    df <- df[order(match(colnames(df), 'ID'))]
    if( !is.null(bind_data) ){
      df <- cbind(df, bind_data)
    }
    #rownames(df) <- x
    
    # add extra ressource
    if( missing(extras) || is_NA(extras) ) return(df)
    extras <- tolower(extras)
    # check which must be links
    with_link <- grepl("^~", extras)
    as_names <- grepl("\\+$", extras)
    as_is <- grepl("^\\.", extras)
    extras <- gsub("^[~.]", "\\1", extras)
    extras <- gsub("\\+$", "", extras)
    extras <- match.arg(extras, several.ok = TRUE)
    
    if( length(i <- which(extras %in% c('ncbi', 'biogps'))) ){
        tolink <- c('ENTREZID', 'Symbol')
        olink <- extras[i[1L]]
        base <- switch(olink
                , biogps = "http://biogps.org/gene/"
                , ncbi = "http://www.ncbi.nlm.nih.gov/gene/")
        .link <- function(x, id = x){
            ok <- !is.na(id) & nzchar(x)
            x[ok] <- sprintf('<a target="_%s" href="%s%s">%s</a>', olink, base, id[ok], x[ok])
            x
        }
        df[tolink] <- sapply(df[tolink], .link, id = ez, simplify = FALSE)
    }
    
    .format <- function(x, i, .link = NULL, id = x, collapse = FALSE){
        
        ok <- !is.na(x)
        sep <- "|" 
        if( do_collapse <- !as_is[i] ){
            is_num <- all(sapply(x[ok], function(...) all(grepl(...)), pattern = '^[0-9]+$'))
            sep <- ifelse(is_num, ",", "\n")
        }
        
        if( with_link[i] && !is.null(.link) ){
            sep <- "<br />\n"
            x[ok] <- mapply(.link, x[ok], id[ok], SIMPLIFY = FALSE)
        } 
        if( do_collapse ) x[ok] <- sapply(x[ok], paste0, collapse = sep)
        if( collapse ) x <- paste0(x, collapse = sep)  
        x
    } 
    
    if( length(i <- which(extras == 'pathway')) ){
        library(reactome.db)
        p <- pid <- bimap_lookup(ez, reactomeEXTID2PATHID)
        # http://pid.nci.nih.gov/search/pathway_landing.shtml?what=graphic&jpg=on&pathway_id=xxxxxx&source=ReactomeImported
        .link <- function(x, id) sprintf('<a target="_pathway" href="http://reactome.org/pathway/%s">%s</a>', id, x)
        if( as_names[i] ){
            p <- sapply(pid, function(x){
                        if( is_NA(x) ) NA
                        else{
                            pn <- bimap_lookup(x, reactomePATHID2NAME)
                            fpn <- .format(pn, i, .link = .link, id = x, collapse = TRUE)
                            df$Pathway <- unlist(fpn)
                        }
                    }, simplify = FALSE)
            .link <- NULL
        }
        
        # format
        df$Pathway <- as.character(.format(p, i, id = pid, .link = .link))
    }
    
    if( length(i <- which(extras == 'kegg')) ){
        k <- bimap_lookup(x, biocann_object('PATH', annotation), multiple = TRUE)
        df$KEGG <- .format(k, i, function(x, id) sprintf('<a target="_pathway" href="http://kegg.org/pathway/%s">%s</a>', x, x))
    }
    
    df
}


#' Checks Gene Symbol Validity
#' 
#' Look for gene symbols that were transformed by Excel nasty default
#' setting, which interpreted them as dates or numeric.
#' 
#' This issue, although well known in the bioinformatics community, regularly 
#' surfaces to the detriment of the person who analyses the data. 
#' It has been reported to affect several publications and even large genomic projects 
#' (\cite{Zeeberg2004;Ziemann2016}).
#' 
#' @param x Character vector of gene symbols to check
#' @param value Logical that indicates if the value of invalid symbols 
#' should be returned (\code{TRUE}), or only their index (\code{FALSE}).
#' @param replace logical that indicates if fixable symbols should be replaced
#' with their correct value.
#' If \code{value = FALSE}, then the corrected symbols are set as names of the 
#' result integer vector, and the incorrect (input) symbols as values.
#' If \code{value = FALSE}, then the corrected symbols are returned as values of the 
#' result character vector, and the incorrect (input) symbols as names. 
#' 
#' @return character or integer vector depending on the value of argument 
#' \var{value}.
#' See help for argument \var{replace}, as this affects what goes into
#' names and values of the result vector.
#' 
#' @export
#' @examples 
#' 
#' bad_symbols <- c('AAA', 'BBB', '1-Dec', '10-JAN')
#' # index of incorrect symbols
#' checkSYMBOLS(bad_symbols)
#' # incorrect symbols
#' checkSYMBOLS(bad_symbols, value = TRUE)
#' # return fixed symbols as names (index vector)
#' checkSYMBOLS(bad_symbols, replace = TRUE)
#' # return fixed symbols as values (mapping)
#' checkSYMBOLS(bad_symbols, value = TRUE, replace = TRUE)
#' 
checkSYMBOLS <- function(x, value = FALSE, replace = FALSE, quiet = TRUE){
	
	# build regexp
	m <- paste0('(', tolower(month.abb), ')', collapse = '|')
	reg <- sprintf("^(([0-9]+)-(%s))|([0-9]\\.[0-9][0-9]E\\+[0-9][0-9])$", m)
	# match
	m <- str_match(x, regex(reg, ignore_case = TRUE))
	# identify matches
	i <- which(!is.na(m[, 1L]))
	if( !length(i) ) return(i)
	
	# distinguish month from scientific numbers
	is_month <- !is.na(m[i, 2L])
	
	if( !quiet ){
		msg <- sprintf("Gene symbols encoded as months: %s\n  Gene symbols encoded as numerics: %s"
						, str_out(m[i[is_month], 1L], total = TRUE)
						, str_out(m[i[!is_month], 1L], total = TRUE)
				)
		warning(msg)
	}
	
	if( !replace ) names(i) <- m[i, 1L] 
	else names(i) <- ifelse(is_month, paste0(m[i, 4L], m[i, 3L]), NA)
	# use value if requested
	if( value ){
		i[] <- m[i, 1L]
		if( replace ) i <- setNames(names(i), i)
	}
	
	# return
	i
}

