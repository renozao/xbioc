% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idtype.R
\name{is_Affymetrix}
\alias{is_Affymetrix}
\alias{is_Illumina}
\alias{is_Agilent}
\title{Testing Identifier Types}
\usage{
is_Affymetrix(x)

is_Illumina(x)

is_Agilent(x)
}
\arguments{
\item{x}{object from which the identifier type is computed by \code{\link{idtype}}.}
}
\description{
This functions tests if identifiers associated with an object is of from a given
type (mainly based on platform).
}
