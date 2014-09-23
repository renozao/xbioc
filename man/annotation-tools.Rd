\name{is.annpkg}
\alias{biocann_mapname}
\alias{biocann_pkgname}
\alias{biocann_pkgobject}
\alias{is.anndb}
\alias{is.annpkg}
\alias{is.orgpkg}
\title{Annotation Tools}
\usage{
  is.annpkg(x)

  is.orgpkg(x)

  is.anndb(x)

  biocann_mapname(annotation, map, all = FALSE)

  biocann_pkgname(x, noext = FALSE)

  biocann_pkgobject(x)
}
\arguments{
  \item{x}{an R object, either a character string or an annotation object.}

  \item{annotation}{names of an annotation package, with \dQuote{.db}
suffix or not.}

  \item{map}{name of a map, e.g., \dQuote{ENTREZID}.}

  \item{all}{logical that indicates if all possible names should be
returned, and only the simple concatenation of the annotation
package's name without \dQuote{.db} and the map name.}

  \item{noext}{logical that indicates if returned package names should
contain the extension '.db'.}
}
\description{
  The \pkg{xbioc} package contains a few utility functions
  to facilitate working with Bioconductor annotations,
  which extends or enhance functions available in packages
  such as \pkg{annotate}.

  \code{is.annpkg} tells if an object is the name of an
  annotation package.

  \code{is.anndb} tells if an object is an annotation db
  object such as \code{hgu133a.db}.

  \code{biocann_mapname} returns the name of a map in an
  annotation package.

  \code{biocann_pkgname} returns the name of an annotation
  package, formated from character strings or extracted
  from annotation objects.

  \code{biocann_pkgobject} retrieve the \code{AnnotationDb}
  object for an annotation package given as a string. The
  package does not need to be already loaded but needs to
  be installed in a library that is in the search path.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

#----------
# is.annpkg
#----------
# check annotation pkg name
is.annpkg('aaa.db')
is.annpkg(c('bbb.db', 'ccc.db'))
is.annpkg(c('ddd', 'eee.db'))

#----------
# is.anndb
#----------
# check AnnotationDb object
library(hgu133plus2.db)
is.anndb(hgu133plus2.db)

#----------
# biocann_mapname
#----------
# build annotation bimap object names
biocann_mapname('hgu133plus2.db', 'ENTREZID')
# '.db' extension is skipped
biocann_mapname('hgu133plus2', 'ENTREZID')
# get all possible map names
biocann_mapname('hgu133plus2.db', 'ENTREZID', all=TRUE)

#----------
# biocann_pkgname
#----------
# annotation package name
biocann_pkgname('aa')
# extract the package name from annotation objects
biocann_pkgname(hgu133plus2.db)
biocann_pkgname(hgu133plus2ENTREZID)

#----------
# biocann_pkgobject
#----------
# get AnnotationDb object
biocann_pkgobject('hgu133plus2') # extension can be skipped
# the package needs to be installed
try( biocann_pkgobject('aaa') )
}

