% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operations.R
\name{log,ExpressionSet-method}
\alias{log,ExpressionSet-method}
\alias{expb}
\alias{expb,numeric-method}
\alias{expb,matrix-method}
\alias{expb,ExpressionSet-method}
\alias{exp,ExpressionSet-method}
\alias{range,ExpressionSet-method}
\alias{quantile.ExpressionSet}
\title{Numeric Computations on ExpressionSet objects}
\usage{
\S4method{log}{ExpressionSet}(x, ...)

expb(x, ...)

\S4method{expb}{numeric}(x, base = exp(1))

\S4method{expb}{matrix}(x, base = exp(1))

\S4method{expb}{ExpressionSet}(x, ...)

\S4method{exp}{ExpressionSet}(x)

\S4method{range}{ExpressionSet}(x, ..., na.rm = FALSE)

\method{quantile}{ExpressionSet}(x, ...)
}
\arguments{
\item{x}{an \code{ExpressionSet} object.}

\item{...}{extra arguments passed to subsequent calls, usually of the corresponding
method in the \pkg{stats} package.}

\item{base}{log base to use.}

\item{na.rm}{logical that indicates if missing values should be omitted from the
computation.}
}
\description{
The \pkg{xbioc} package defines some generics and methods to apply numeric transformations
to \code{ExpressionSet} objects, which is convenient when working on gene expression
deconvolution algorithms, where scale (log/linear) may matter.

\code{log} log-transforms the expression matrix of \code{\link{ExpressionSet}} objects.
}
