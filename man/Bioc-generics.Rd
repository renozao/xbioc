\docType{methods}
\name{featureNames}
\alias{exprs}
\alias{exprs,matrix-method}
\alias{exprs-methods}
\alias{featureNames}
\alias{featureNames<-}
\alias{featureNames<-,matrix-method}
\alias{featureNames,matrix-method}
\alias{featureNames<--methods}
\alias{featureNames-methods}
\alias{sampleNames}
\alias{sampleNames<-}
\alias{sampleNames,matrix,ANY}
\alias{sampleNames<-,matrix-method}
\alias{sampleNames,matrix-method}
\alias{sampleNames<--methods}
\alias{sampleNames-methods}
\title{Extracting Feature Names}
\usage{
  featureNames(object)

  featureNames(object)<-value

  sampleNames(object)

  sampleNames(object)<-value

  exprs(object)

  \S4method{exprs}{matrix}(object)
}
\arguments{
  \item{object}{Object, possibly derived from class \code{eSet}.}

  \item{value}{Character vector containing feature or sample names.}
}
\description{
  The \pkg{xbioc} package provides extra methods for the
  generic \code{\link[Biobase]{featureNames}} and
  \code{\link[Biobase]{sampleNames}}, that extend the
  original Bioconductor interface.
}
\section{Methods}{
\describe{

\item{exprs}{\code{signature(object = "matrix")}: Simply returns \code{object}.
This method is defined so that function can seamlessly handle both matrix and
\code{ExpressionSet} objects.
}


\item{featureNames}{\code{signature(object = "matrix")}: Returns the row names of \code{object}.
}


\item{featureNames<-}{\code{signature(object = "matrix")}: Set the row names of \code{object}.
}


\item{sampleNames}{\code{signature(object = "matrix")}: Returns the column names of \code{object}.
}


\item{sampleNames<-}{\code{signature(object = "matrix")}: Set the row names of \code{object}.
}

}

}
\keyword{methods}

