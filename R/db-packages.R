# Project: xbioc
# 
# Function for annotation packages
#
# Author: Renaud Gaujoux
# Created: Jul 3, 2014
###############################################################################

#' Annotation Tools
#' 
#' @description
#' The \pkg{CellMix} package contains a few utility functions to facilitate 
#' working with Bioconductor annotations, which extends or enhance functions
#' available in packages such as \pkg{annotate}.
#' 
#' \code{is.annpkg} tells if an object is the name of an annotation package.
#' 
#' @param x an R object, either a character string or an annotation object.
#' 
#' @rdname annotation-tools
#' @export
#' @examples
#' # check annotation pkg name
#' is.annpkg('aaa.db')
#' is.annpkg(c('bbb.db', 'ccc.db'))
#' is.annpkg(c('ddd', 'eee.db'))
is.annpkg <- function(x) is.character(x) && length(x)>0L && all(grepl("\\.db$", x))

#' \code{is.anndb} tells if an object is an annotation db object such as 
#' \code{hgu133a.db}.
#' 
#' @rdname annotation-tools
#' @export
#' @examples
#' # check AnnotationDb object
#' library(hgu133plus2.db)
#' is.anndb(hgu133plus2.db)
is.anndb <- function(x) is(x, 'AnnotationDb')

#' \code{biocann_mapname} returns the name of a map in an annotation package.
#' 
#' @param annotation names of an annotation package, with \dQuote{.db} 
#' suffix or not.
#' @param map name of a map, e.g., \dQuote{ENTREZID}.
#' @param all logical that indicates if all possible names should be 
#' returned, and only the simple concatenation of the annotation 
#' package's name without \dQuote{.db} and the map name.
#' 
#' @rdname annotation-tools
#' @export
#' @examples
#' # build annotation bimap object names
#' biocann_mapname('hgu133plus2.db', 'ENTREZID')
#' # '.db' extension is skipped
#' biocann_mapname('hgu133plus2', 'ENTREZID')
#' # get all possible map names
#' biocann_mapname('hgu133plus2.db', 'ENTREZID', all=TRUE)
biocann_mapname <- function(annotation, map, all=FALSE){
    
    base <- biocann_pkgname(annotation, noext=TRUE)
    sep <- ''
    if( all ) sep <- c(sep, '2')
    paste(base, sep, map, sep='')
    
}
#' \code{biocann_pkgname} returns the name of an annotation package, formated from character strings
#' or extracted from annotation objects.
#' 
#' @param noext logical that indicates if returned package names should 
#' contain the extension '.db'.
#' 
#' @rdname annotation-tools
#' @export
#' @examples
#' # annotation package name
#' biocann_pkgname('aa')
#' # extract the package name from annotation objects
#' biocann_pkgname(hgu133plus2.db)
#' biocann_pkgname(hgu133plus2ENTREZID)
biocann_pkgname <- function(x, noext=FALSE){
    
    # extract name from annotation db object
    if( is.anndb(x) ) x <- x$packageName
    if( is(x, 'ProbeAnnDbBimap' ) ) x <- strsplit(x@objTarget, ' ')[[1]][2]
    
    if( !is.character(x) )
        stop("Invalid argument `x`: character string expected [", class(x), ']')
    
    base <- sub("\\.db$", "", x)
    if( noext ) base 
    else ifelse(nchar(base), paste(base, ".db", sep=''), '')
    
}
#' \code{biocann_pkgobject} retrieve the \code{AnnotationDb} object for an annotation
#' package given as a string.
#' The package does not need to be already loaded but needs to be installed in a library that 
#' is in the search path.
#' 
#' @rdname annotation-tools
#' @export
#' @examples
#' 
#' # get AnnotationDb object
#' biocann_pkgobject('hgu133plus2') # extension can be skipped
#' # the package needs to be installed
#' try( biocann_pkgobject('aaa') )
#' 
biocann_pkgobject <- function(x){
    
    pkg <- biocann_pkgname(x)
    get(pkg, envir=asNamespace(pkg))
    
}

#' Retrieving Bioconductor Annotation Maps
#' 
#' The function \code{biocann_object} retrieves annotation objects, like bimaps, from 
#' Bioconductor annotation packages.
#' It is similar to the function \code{getAnnMap} from the \pkg{annotate} package, 
#' except that it also accepts annotation -- bimap -- objects, 
#' and will try to install missing packages if not found (see section Details).
#' 
#' If an annotation package is specified as a character string, and is not found in the 
#' search path, and if R runs in interactive mode, then the user is asked whether one 
#' should try install the missing package.
#' Default response is 'yes' so that simply hitting return will install the package
#' via \code{\link{install.packages}} and load it.
#' An error is thrown if this eventually fails.
#' 
#' @param to target annotation field as a character string, e.g., \dQuote{ENTREZID}, 
#' \dQuote{ENSEMBL}, or an annotation package or db which means that one wants to 
#' retrieve a mapping to its corresponding primary identifier. 
#' If \code{from} is missing, \code{to} must be the name of an annotation package, 
#' i.e. ends with \dQuote{db}), in which case it tries loading the package and return 
#' the whole annotation db object; or any annotation package db or map object such as 
#' \code{AnnDbBimap}, \code{ChipDb} or \code{OrgDb} objects, which are returned unchanged.
#' @param from source annotation package as a character string e.g. \code{"hgu133a.db"}.
#' @param optional logical that indicates if the function should return \code{NULL} if the 
#' mapping cannot be found (\code{TRUE}), or throw an error.
#' Note that this does not apply to the installation part: if a required annotation 
#' package is missing, an error is thrown even if \code{optional=TRUE}. 
#' 
#' @return a \code{\link{ProbeAnnDbBimap}} if annotation is not missing, 
#' a \code{ProbeAnnDb} object otherwise.
#' 
#' @export
#' @examples
#' 
#' # db package object
#' biocann_object('hgu133plus2.db')
#' 
#' # bimap from hgu133plus2 probe id to ENTREZID 
#' biocann_object('ENTREZID', 'hgu133plus2.db')
#' 
#' # reversed bimap from UNIGENE to hgu133plus2 probe id
#' biocann_object('hgu133plus2.db', 'UNIGENE')
#' # this is equivalent to using the annotation package object (no quotes),
#' # when the package is already loaded (=> helpful in interactive session with auto-completion)
#' biocann_object(hgu133plus2.db, 'UNIGENE')
#' 
biocann_object <- function(to, from=NULL, optional=FALSE){
    
    if( is.null(from) ){
        # simply return if the argument is already an annotation Bimap
        if( is(to, 'AnnDbBimap') 
                || is(to, 'ChipDb') 
                || is(to, 'OrgDb') ){
            return(to)
        }
        # to must be the name of an annotation package 
        if( !is.annpkg(to) )
            stop("Invalid annotation package name ['",to,"']")
        
        # try load and ask if not possible 
        if( !uq_requirePackage(to, load=TRUE) )
            stop("Aborted computation due to missing annotation package '", to, "'.")
        
        ns <- asNamespace(to)
        return(get(to, ns))
        
    }
    
    # use annotation string from ExpresssionSet objects
    if( isExpressionSet(from) ) from <- annotation(from)
    
    if( is.character(from) ){
        
        if( is.annpkg(to) || is.anndb(to) ){ # e.g., we want 'ENTREZID -> pkg.db' 
            if( is.annpkg(from) )
                stop("Can only map between an annotation package and a data field.",
                        " [map='", to, "' - annotation=", class(from), "]")
            
            # get map and reverse it
            m <- biocann_object(from, to, optional = optional)
            if( !is.null(m) ) return( revmap(m) )
            else return()
        }
        
        # load db package
        annenv <- biocann_object(from)
        
    }else if( is.environment(from) ){
        annenv <- from
    }
    
    if( !is.character(to) )
        stop("Invalid argument `to`: expected character string [", class(to), '] when `from` is [', class(from), ']')
    
    # get all potential map names 
    maps <- biocann_mapname(from, to, all=TRUE)
    for( mname in maps ){
        if( exists(mname, annenv) )	return( get(mname, annenv) )
    }
    
    # error if not optional
    if( !optional ){
        annpkg <- biocann_pkgname(from)
        stop("Could not find map for '", to, "' in package '", annpkg, "'")
    }
    NULL
}

####################
# ORGANISM PACKAGES
####################
#' Bioconductor Organism Data and Packages
#' 
#' @param organism Organism name (case insensitive).
#' Can be specified using the common name (e.g., rat), a partial match of the latin name (e.g., Mus [musculus]), 
#' or its abbreviation (e.g., Mm). 
#' 
#' The following organisms are currently supported:
#' 
#' human, arabidopsis, bovine, canine, chicken, chimp, malaria, mouse, pig, rat, rhesus, worm, zebrafish
#'  
#' @export
#' @examples
#'
#' # list all 
#' biocann_orgdb()
#' 
#' # common name
#' biocann_orgdb('human')
#' 
#' # latin name
#' biocann_orgdb('canis')
#' # any partial match works
#' biocann_orgdb('canis lupus')
#' 
#' # abbreviation
#' biocann_orgdb('Mm')
#' 
biocann_orgdb <- local({
    
    x <- list(human = 'Homo sapiens'
        , mouse = 'Mus musculus' 
        , rat = 'Rattus norvegicus'
        , arabidopsis = 'Arabidopsis thaliana'
        , chicken = 'Gallus gallus'
        , canine = 'Canis lupus familiaris'
        , bovine = 'Bos taurus'
        , rhesus = 'Macaca mulatta'
        , chimp = 'Pan troglodytes'
        , pig = 'Sus scrofa'
        , worm = 'Caenorhabditis elegans'
        , zebrafish = 'Danio rerio'
        , malaria = 'Plasmodium falciparum')
    x <- data.frame(organism = t(data.frame(x)), stringsAsFactors = FALSE)
    
    # add abbreviations
    x$abbrv <- sapply(strsplit(x$organism, ' '), function(x) paste0(substr(x,1,1), collapse = ''))
    x[c('canine', 'rhesus'), 'abbrv'] <- c('Cf', 'Mmu')
    
    # db0 package name
    x$org.db0 <- sprintf('%s0', paste0(rownames(x), '.db'))
    # org package name
    x$org.db <- sprintf('org.%s.eg.db', x$abbrv)
    spe_org <- c(arabidopsis='tair', malaria='plasmo')
    x[names(spe_org), 'org.db'] <- sapply(seq_along(spe_org), function(i) gsub("\\.eg\\.", paste0(".", spe_org[i], "."), x[names(spe_org)[i], 'org.db'])) 
    
    # reorder
    x$db0 <- rownames(x)
    .map <- x[order(rownames(x)), ]
    
    function(organism){
        if( missing(organism) ) return(.map)
        x <- organism
        
        if( x %in% rownames(.map) ) i <- x
        else if( is.na(i <- match(toupper(x), toupper(.map$organism))) ){
                if( !length(i <- grep(paste0("^", x), .map$organism, ignore.case = TRUE)) ){
                        if( is.na(i <- match(toupper(x), toupper(.map$abbrv))) ){
                                stop("Invalid organism: ", x, "\n  Must be one of: ", str_out(paste0(rownames(.map), " [", .map$abbrv, ']'), Inf, quote = FALSE))
                        }
                }
        }
        .map[i, , drop = FALSE]
    }
})




#' Looking Up Keys in Identifier Maps
#' 
#' Look up for keys in a map.
#' 
#' @param keys keys to look up
#' @param map Map as a bimap object, e.g., \code{hgu133plus2ENTREZID}.
#' @param multiple logical that indicates if all matches should be included, 
#' or if the result list should be reduced to a one-to-one mapping vector.
#' Possible values are:
#' \itemize{
#' \item \code{TRUE} result is a list with potential multiple matches
#' \item \code{NA} results is a vector with multiple matches set to NA 
#' \item \code{FALSE} or \code{'first'} results is a vector with multiple matches resolved by choosing
#' the \emph{first} match.
#' \item \code{'last'} results is a vector with multiple matches resolved by choosing
#' the \emph{last} match.
#' }
#' 
#' @export
bimap_lookup <- function(keys, map, multiple = TRUE){
    res <- if( !is.list(map) ){
        AnnotationDbi::mget(keys, map, ifnotfound=NA)
    }else{
        res <- NAmap(keys) #setNames(as.list(rep(NA, length(keys))), keys)
        # early exit on empty map
        if( !length(map) ) return(res)
        mk <- keys[keys %in% mappedkeys(map)]
        res[mk] <- map[mk]
        res
    }

    # process multiple matches
    if( !isTRUE(multiple) ){
        if( is_NA(multiple) ){ # set multiple match to NA
            res[sapply(res, length) > 1L] <- NA
            res <- unlist(res)
        }else if( isFALSE(multiple) || multiple == 'first' ) # choose first 
            res <- sapply(res, head, 1L)
        else if( multiple == 'last' ) # choose last 
            res <- sapply(res, tail, 1L)
        else stop("Invalid argument `multiple`: must be TRUE, FALSE, NA, first' or 'last'")
    }
    
    # return
    res
}
