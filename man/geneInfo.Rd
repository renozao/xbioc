\name{geneInfo}
\alias{geneInfo}
\title{Simple Feature Annotation}
\usage{
  geneInfo(x, organism = "human")
}
\arguments{
  \item{x}{ENTREZ gene ids}

  \item{organism}{Organism name (case insensitive). Can be
  specified using the common name (e.g., mouse), a partial
  match of the latin name (e.g., Mus [musculus]), or its
  abbreviation (e.g., Mm).

  The following organisms are currently supported:

  human, arabidopsis, bovine, canine, chicken, chimp,
  malaria, mouse, pig, rat, rhesus, worm, zebrafish}
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

