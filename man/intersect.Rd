\docType{methods}
\name{intersect}
\alias{intersect}
\alias{intersect,MatrixData,character-method}
\alias{intersect,MatrixData,logical-method}
\alias{intersect,MatrixData,MatrixData-method}
\title{Enhanced Subsetting for Matrix-like Data}
\usage{
  \S4method{intersect}{MatrixData,logical}(x, y)

  \S4method{intersect}{MatrixData,character}(x, y)

  \S4method{intersect}{MatrixData,MatrixData}(x, y)
}
\arguments{
  \item{x}{vectors (of the same mode) containing a sequence
    of items (conceptually) with no duplicated values.}

  \item{y}{vectors (of the same mode) containing a sequence
    of items (conceptually) with no duplicated values.}
}
\description{
  These methods subset Matrix-like data only keeping rows
  whose names are in common with some other data.
}
\section{Methods}{
\describe{

\item{intersect}{\code{signature(x = "MatrixData", y = "logical")}: This method is equivalent to \code{\link[base]{subset}(x, y)}.
}

\item{intersect}{\code{signature(x = "MatrixData", y = "character")}: Subset a matrix-like object by only keeping the rows whose
names are in a given reference character vector.
}

\item{intersect}{\code{signature(x = "MatrixData", y = "MatrixData")}: Subset a matrix-like object by only keeping the rows whose
names are in common with another matrix-like data.

This is a shortcut for \code{intersect(x, featureNames(y), ...)}.
}

}

}
\keyword{methods}

