# Project: xbioc
# 
# Author: Renaud Gaujoux
# Created: Jul 20, 2014
###############################################################################

#' @include AllGenerics.R
NULL

#' Identifying Gene or Probe ID Type
#' 
#' The S4 generic \code{idtype} automatically determines the type 
#' of gene/feature identifiers stored in objects, based on a combination of 
#' regular expression patterns and test functions.
#' 
#' It uses a heuristic based on a set of regular expressions and functions 
#' that uniquely match most common types of identifiers, such as Unigene, 
#' Entrez gene, Affymetrix probe ids, Illumina probe ids, etc..
#' 
#' @param object an R object that contains the gene identifiers whose type is to 
#' be determined.
#' @param ... extra argument to allow extension, generally passed down to 
#' \code{idtype,character-method}.
#' See each method's description for more details.
#' 
#' @export
#' @inline
setGeneric('idtype', function(object, ...) standardGeneric('idtype'))

#' \code{is.idtype} tells if a given character vector contains valid types.
#' 
#' @rdname biocIDs
#' @export
is.idtype <- function(type){
    if( !is.character(type) ) FALSE
    else type %in% idtype()
}

is.multiple <- function(x){
    sapply(x, length) > 1L
}

#' Utility function for Biological Identifiers
#' 
#' \code{is.probeid} tells if given IDs are probe IDs.
#' 
#' @param x an R object
#' @param ... extra arguments passed to \code{\link{idtype}}.
#' 
#' @family biocIDs
#' @rdname biocIDs
#' @export
is.probeid <- function(x, ...) is.probetype(idtype(x, ...))
#' \code{is.probetype} tells if given types are probe ID types.
#' 
#' @param type an identifier type as returned by \code{\link{idtype}}.
#' 
#' @rdname biocIDs
#' @export
is.probetype <- function(type) type == 'ANNOTATION' | (is.idtype(type) & grepl("^\\.", type))

#' Method for when \code{idtype} is called with its first argument missing, 
#' in which case it returns all or a subset of the known type names as a character 
#' vector, or optionally as a list that contains their definition, i.e. a regular 
#' expression or a matching function.
#' 
#' @param def a logical or a subsetting vector, used when \code{object} is missing, 
#' which indicates that the result should contain the definition of the matching 
#' pattern/function of each type, or which type's deifnition should be included
#' in the result list.
#'   
#' @examples 
#' 
#' # all known types
#' idtype()
#' # with their definitions
#' idtype(def=TRUE)
#' idtype(def='ENTREZID')
#' idtype(def=c('ENTREZID', 'ENSEMBLTRANS'))
#' 
setMethod('idtype', 'missing' 
        , local({
                    
        # type matching patterns/functions
        .defs <- list(
                UNIGENE="^[A-Z][a-z]\\.[0-9]+$"
                , ENSEMBL= NA_character_
                , ENSEMBLTRANS= NA_character_
                , ENSEMBLPROT= NA_character_
                , ENTREZID="^[0-9]+$"
                , IMAGE = "^IMAGE:[0-9]+$"
                , GOID="^GO:[0-9]+$"
                , PFAM="^PF[0-9]+$"
                , REFSEQ="^[XYN][MPR]_[0-9.]+$"
                , ENZYME="^[0-9]+(\\.(([0-9]+)|-)+){3}$"
                , MAP="^(([0-9]{1,2})|([XY]))((([pq])|(cen))(([0-9]+(\\.[0-9]+)?)|(ter))?(-([0-9]{1,2})|([XY]))?(([pq]?)|(cen))((ter)|([0-9]+(\\.[0-9]+)?))?)?)?$"
                , GENEBANK=c("^[A-Z][0-9]{5}$", "^[A-Z]{2}[0-9]{6}$")
                , GENEBANK="^[A-Z]{3}[0-9]{5}$"
                , GENEBANK="^[A-Z]{4}[0-9]{8}[0-9]?[0-9]?$"
                , GENEBANK="^[A-Z]{5}[0-9]{7}$"
                , GENENAME = " "
                , .Affymetrix="(^AFFX[-_])|(^[0-9]+_([abfgilrsx]_)?([as]t)|(i))$"
                , .Illumina="^ILMN_[0-9]+$"
                , .Agilent = "^A_[0-9]+_P[0-9]+$"
                , .nuID=function(x) !is.na(nuIDdecode(x, error=NA))
        #, SYMBOL="^[-0-9a-zA-Z':@.]*[a-zA-Z][-0-9a-zA-Z':@.]*$"
        )
        
        # actual function
        function(object, def=FALSE){
            # load ENSEMBL prefixes from data table if needed
            if( is_NA(.defs$ENSEMBL) ){
              d <- ldata('ENSEMBL_prefix', package = 'xbioc')
              pref <- paste0(as.character(d$Prefix), collapse = '|')
              .defs$ENSEMBL <<- sprintf("^(%s)G[0-9]+(\\.[0-9]+.*)?$", pref)
              .defs$ENSEMBLTRANS <<- sprintf("^(%s)T[0-9]+(\\.[0-9]+.*)?$", pref)
              .defs$ENSEMBLPROT <<- sprintf("^(%s)P[0-9]+(\\.[0-9]+.*)?$", pref)
              
            }
            
            if( isFALSE(def) ) unique(names(.defs))
            else if( isTRUE(def) ) .defs
            else if( is.vector(def) ){
                if( is.character(def) && length(def) == 1L ){
                    res <- .defs[names(.defs) %in% def]
                    if( length(res) == 1L ) res[[1L]]
                    else res
                }else .defs[def]
            }else
                stop("idtype - Invalid argument `def` of class:'",class(def),"': expecting FALSE, TRUE or a subsetting vector.")
        }
    })
)
#' Detects the type of identifiers used in the row names of a matrix. 
setMethod('idtype', 'matrix', function(object, ...) idtype(rownames(object), ...) )
#' Dummy method -- defined for convenience -- that returns \code{''}
setMethod('idtype', 'NULL', function(object, ...) '' )

#' This is the workhorse method that determine the type of ids contained in 
#' a character vector. 
#' 
#' @param each logical indicating whether the type of each element should be 
#' returned (\code{TRUE}) or only the type of the vector as a whole (default).
#' @param limit specification for limiting which elements are used to 
#' detect the type of identifiers.
#' If a single numeric, then only the first \code{limit} elements 
#' are used. Otherwise it must be a subsetting logical or numeric vector.
#' @param no.match character string that specifies the string to use when the 
#' type cannot be determined.    
#' 
#' The IDs can be either:
#' \itemize{
#' \item{probe IDs (e.g. 123456_at or ILMN_123456 for Affymetrix or Illumina 
#' chips respectively), the type starts with a dot \code{'.'}, allowing the 
#' subsequent handling of such IDs as a group.}
#' \item{other biological ID types, the result are character strings such as 
#' those used as attributes in Bioconductor annotation packages 
#' (e.g. \code{"ENTREZID"} or \code{"ENSEMBL"})}
#' \item{Names of annotation packages e.g. \code{"hgu133plus2.db"}.}
#' }
#' 
#' This function is able to identify the following ID types using regular 
#' expression patterns or dedicated function:
#' \itemize{
#' \item ENSEMBL = "^ENSG[0-9]+$"
#' \item ENSEMBLTRANS = "^ENST[0-9]+$"
#' \item ENSEMBLPROT = "^ENSP[0-9]+$"
#' \item ENTREZID = "^[0-9]+$"
#' \item IMAGE = "^IMAGE:[0-9]+$"
#' \item GOID = "^GO:[0-9]+$"
#' \item PFAM = "^PF[0-9]+$"
#' \item REFSEQ = "^N[MP]_[0-9]+$"
#' \item ENZYME = "^[0-9]+(\\.(([0-9]+)|-)+){3}$"
#' \item MAP = "^[0-9XY]+((([pq])|(cen))(([0-9]+(\\.[0-9]+)?)|(ter))?(-([0-9XY]+)?(([pq]?)|(cen))((ter)|([0-9]+(\\.[0-9]+)?))?)?)?$"
#' \item GENEBANK (Nucleotide) = "^[A-Z][0-9]{5}$" | "^[A-Z]{2}[0-9]{6}$" 
#' \item GENEBANK (Protein) = "^[A-Z]{3}[0-9]{5}$"
#' \item GENEBANK (WGS) = "^[A-Z]{4}[0-9]{8}[0-9]?[0-9]?$"
#' \item GENEBANK (MGA) = "^[A-Z]{5}[0-9]{7}$"
#' \item GENENAME = " "
#' \item .Affymetrix = "(^AFFX-)|(^[0-9]+_([abfgilrsx]_)?([as]t)|(i))$"
#' \item .Illumina = "^ILMN_[0-9]+$"
#' \item .Agilent = "^A_[0-9]+_P[0-9]+$"
#' \item .nuID = use the function \code{\link{nuIDdecode}} to try converting the 
#' ids into nucleotide sequences. Identification is positive if no error is 
#' thrown during the conversion.  
#' }
#' 
#' @return a single character string (possibly empty) if \code{each=FALSE} (default) 
#' or a character vector of the same "length" as \code{object} otherwise.
#' @noMd
#' @export
#' @examples
#' 
#' idtype("12345_at")
#' idtype(c("12345_at", "23232_at", "555_x_at"))
#' # mixed types
#' ids <- c("12345_at", "23232_at", "Hs.1213")
#' idtype(ids) # not detected
#' idtype(ids, each=TRUE)
#' 
setMethod('idtype', 'vector',
        function(object, each=FALSE, limit=NULL, no.match=''){
            
            no.match <- as.character(no.match)
            
            if( length(object) == 0 ) return(no.match)
            # handle annotation package names
            if( is.character(object) && all(is.annpkg(object)) ){
                res <- 
                        if( length(object) == 1L ) idtype(biocann_object(object))
                        else sapply(object, function(x){ idtype(biocann_object(x), no.match=NA) })
                return(res)
            }
            
            if( each ) res_vec <- setNames(rep(no.match, length(object)), as.character(object)) 
            
            ids <- object
            # convert to character if necessary
            if( !is.numeric(object) ) ids <- as.character(ids)
            # remove missing
            i_process <- which(!is.na(ids))
            # remove empty strings
            if( is.character(ids) ) i_process <- i_process[nzchar(ids[i_process])]
            
            # early exit if empty
            if( length(i_process) == 0 ){
                if( each ) return(res_vec)
                else return(no.match)
            }
            
            # use only a subset of ids 
            if( !is.null(limit) ){
                if( !is.numeric(limit) && !is.logical(limit) )
                    stop("idtype - Invalid value for argument `limit`: must be numeric or logical vector.")
                
                i_process <- if( length(limit) == 1L ) head(i_process, limit)
                            else if( is.logical(limit) ) i_process[limit]
                            else i_process[as.integer(limit)]
            }
            
            # subset original ids 
            ids <- ids[i_process]
            if( length(ids) == 0 ){
                warning("idtype - Input vector is empty or filled wih NAs or empty strings: returned NA(s).")
                if( each ) return( res_vec)                    
                else return(as.character(NA))
            }
            patterns <- idtype(def=TRUE)		
            
            # special case for integer vectors
            if( is.integer(ids) ){
                if( each ){
                    res_vec[i_process] <- setNames(rep("_INDEX_", length(ids)), ids)
                    return(res_vec)
                }else return("_INDEX_")
            }
            
            .check <- function(def, keys){
                if( is.character(def) ){
                    res <- grepl(def[1L], keys)
                    if( length(def)>1L ){
                        for( i in 2:length(def) ){
                            res <- res | grepl(def[i], keys)
                        }
                    }
                    res
                } else if( is.function(def) ) def(keys)
                else stop("Invalid ID matcher object [", class(def), ']')
            }
            
            if( each ){
                res <- setNames(rep(no.match, length(ids)), ids)
                idx <- seq_along(res)
            }
            for( i in seq_along(patterns) ){
                p <- patterns[[i]]
                
                if( !each ){
                    if( all(.check(p, head(ids, 3))) && all(.check(p, ids)) )
                        return(names(patterns)[i])
                }else{
                    t <- .check(p, ids[idx])
                    ok <- which(t)
                    if( length(ok) > 0 ){
                        res[idx[ok]] <- names(patterns)[i]
                        idx <- idx[-ok]
                        if( length(idx) == 0L ) break
                    }
                }
            }
            
            if( each ){
                res_vec[i_process] <-  res
                res_vec
            }else no.match
        }
)

#' Detects the type of identifiers used in the feature names of an \code{ExpressionSet} object. 
setMethod('idtype', 'ExpressionSet', function(object, ...) idtype(featureNames(object), ...) )
#' Detects the type of the primary identifiers of a probe annotation bimap object.
#'  
#' To speedup the identification, only the first 500 probes are used by default, 
#' since the IDs are very likely to have been curated and to be of the same type.
#' This can be changed using argument \code{limit}.
#' 
setMethod('idtype', 'ProbeAnnDbBimap', function(object, limit=500L, ...) idtype(keys(object), limit=limit, ...) )
#' Detects the type of the identifiers of a chip annotation object.
#' 
#' To speedup the identification, only the first 500 probes are used by default, 
#' since the IDs are very likely to have been curated and to be of the same type.
#' This can be changed using argument \code{limit}.
setMethod('idtype', 'ChipDb', function(object, limit=500L, ...) idtype(keys(object), limit=limit, ...) )
#' Detects the type of the identifiers of an organism annotation object.
#' 
#' To speedup the identification, only the first 500 probes are used by default, 
#' since the IDs are very likely to have been curated and to be of the same type.
#' This can be changed using argument \code{limit}.
setMethod('idtype', 'AnnDbBimap', function(object, limit=500L, ...) idtype(keys(object), limit=limit, ...) )
#' Detects the type of all elements in a list, but provides the option of 
#' detecting the type of each element separately.
#' 
setMethod('idtype', 'list', function(object, ...) idtype(unlist(object), ...) )
setMethod('idtype', 'data.frame', function(object, ...) sapply(object, idtype, ...) )

#' Dummy method -- defined for convenience -- that returns \code{''}
setMethod('idtype', 'NULL', function(object, ...) '' )

#' Testing Identifier Types
#' 
#' This functions tests if identifiers associated with an object is of from a given 
#' type (mainly based on platform).
#' 
#' @param x object from which the identifier type is computed by \code{\link{idtype}}.  
#' 
#' @rdname idtype-test
#' @export
is_Affymetrix <- function(x){ idtype(x) == '.Affymetrix'}
#' @rdname idtype-test
#' @export
is_Illumina <- function(x){ idtype(x) == '.Illumina'}
#' @rdname idtype-test
#' @export
is_Agilent <- function(x){ idtype(x) == '.Agilent'}

#' Convert nuID to Nucleotide Sequence
#' 
#' The function \code{nuIDdecode} converts a nuID string back to 
#' the nucleotide sequence it encodes.
#' nuIDs are identifiers used as primary keys in some Illumina annotation packages, 
#' and are based on a hash of the probe sequence itself. 
#' 
#' This function is an adapted version of the \code{lumi::id2seq} from the 
#' \pkg{lumi} package, so that it can throws errors instead of warnings.
#' It is used in \code{\link{idtype}} to infer the type of nuID vectors.
#' 
#' @param id nuID character vector
#' @param error a logical indicating whether an error should be thrown in case 
#' of invalid input nuID (default) or only a warning, as in the original function
#' \code{lumi::id2seq}.
#' 
#' @source function \code{lumi::id2seq} in the \pkg{lumi} package
#' 
#' @author Pan Du
#' 
#' Adaptation for \pkg{xbioc}: Renaud Gaujoux
#' @cite Du2007
#' 
#' @export
#' @examples
#' nuIDdecode('XERDqYYc2A')
#' try(nuIDdecode('XERDqYYc2F'))
#' nuIDdecode('XERDqYYc2F', error=FALSE)
#' 
nuIDdecode <- function (id, error = TRUE) 
{
    if (length(id) > 1) {
        return(sapply(id, nuIDdecode, error=error))
    }
    else {
        code <- c(LETTERS, letters, as.character(0:9), "_", ".")
        ind <- 1:length(code)
        names(ind) <- code
        if (length(grep("[^a-zA-Z0-9_.]", id)) > 0){			
            msg <- "Input id is not a nuID!"
            if( isTRUE(error) ) stop(msg)
            else return(NA) 
            #else warning(msg)
        }
        id <- substring(id, 1:nchar(id), 1:nchar(id))
        num <- as.numeric(ind[id]) - 1
        checkCode <- num[1]
        num <- num[-1]
        if (checkCode == 63){
            msg <- "Coding error or not a nuID!\n Check code should not include \".\"!"
            if( isTRUE(error) ) stop(msg)
            else if( is.na(error) ) return(NA) 
            else warning(msg)
        }
        cutLen <- checkCode%%3
        res <- floor(checkCode/3)
        codon <- rbind(floor(num/4^2), floor((num%%4^2)/4), num%%4)
        checkSum <- sum(num)
        if (res != checkSum%%21) {
            msg <- "Coding error or not a nuID!"
            if( isTRUE(error) ) stop(msg)
            else if( is.na(error) ) return(NA)
            else warning(msg)
        }
        nucleotide <- c("A", "C", "G", "T")
        codon <- nucleotide[codon + 1]
        len <- length(codon)
        if (len <= cutLen) { # added safeguard
          msg <- "Coding error or not a nuID!"
          if( isTRUE(error) ) stop(msg)
          else if( is.na(error) ) return(NA)
          else warning(msg)
        }
        seq <- paste(codon[1:(len - cutLen)], collapse = "")
        return(seq)
    }
}
