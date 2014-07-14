\name{geneInfo}
\alias{geneInfo}
\title{Simple Feature Annotation}
\usage{
  geneInfo(x, annotation = "human")
}
\arguments{
  \item{x}{ENTREZ gene ids}

  \item{annotation}{Annotation package to use. It can be
  the name of an annotation package or of an organism
  supported by \code{\link{biocann_orgdb}}.}
}
\description{
  Simple Feature Annotation
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

geneInfo(1:20)
}
\seealso{
  \code{\link{biocann_orgdb}}
}

