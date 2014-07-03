\docType{methods}
\name{compareFeatures}
\alias{compareFeatures}
\alias{compareFeatures,ANY,ANY-method}
\alias{compareFeatures,ANY,ChipDb-method}
\alias{compareFeatures,ANY,matrix-method}
\alias{compareFeatures,character,character-method}
\alias{compareFeatures,ChipDb,ANY-method}
\alias{compareFeatures,matrix,ANY-method}
\alias{compareFeatures-methods}
\title{Compare Features Sets}
\usage{
  compareFeatures(x, y, ...)
}
\arguments{
  \item{x,y}{objects with associated featureNames}

  \item{...}{extra arguments eventually passed to
  \code{compareFeatures,character,character}.}
}
\description{
  Compares the set of features associated with an object,
  e.g., \code{\link[Biobase]{ExpressionSet}} or
  \code{\link[AnnotationDbi]{ChipDb-class}} objects.
}
\section{Methods}{
  \describe{

  \item{compareFeatures}{\code{signature(x = "character", y
  = "character")}: Workhorse function comparing the feature
  sets \code{x} and \code{y}. }

  }
}
\keyword{methods}

