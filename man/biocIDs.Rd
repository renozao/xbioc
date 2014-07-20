\name{is.idtype}
\alias{is.idtype}
\alias{is.probeid}
\alias{is.probetype}
\title{\code{is.idtype} tells if a given character vector contains valid types.}
\usage{
  is.idtype(type)

  is.probeid(x, ...)

  is.probetype(type)
}
\arguments{
  \item{x}{an R object}

  \item{...}{extra arguments passed to
  \code{\link{idtype}}.}

  \item{type}{an identifier type as returned by
  \code{\link{idtype}}.}
}
\description{
  \code{is.idtype} tells if a given character vector
  contains valid types.

  \code{is.probeid} tells if given IDs are probe IDs.

  \code{is.probetype} tells if given types are probe ID
  types.
}

