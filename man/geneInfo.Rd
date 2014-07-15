\name{geneInfo}
\alias{geneInfo}
\title{Simple Feature Annotation}
\usage{
  geneInfo(x, annotation = "human",
    link = c("none", "bioGPS", "NCBI"))
}
\arguments{
  \item{x}{ENTREZ gene ids}

  \item{annotation}{Annotation package to use. It can be
  the name of an annotation package or of an organism
  supported by \code{\link{biocann_orgdb}}.}

  \item{link}{Indicates the type of online resource to link
  gene identifiers. If not \code{'none'}, then ENTREZID and
  Symbols columns contain HTML links (<a> tags) that links
  to the corresponding gene page.}
}
\description{
  Simple Feature Annotation
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# Entrez IDs
geneInfo(1:5)
geneInfo(1:5, link = 'bioGPS')

# probe ids
geneInfo(c("1007_s_at", "1053_at", "117_at", "121_at", "not_a_probe_id"), 'hgu133plus2.db')
geneInfo(c("1007_s_at", "1053_at", "117_at", "121_at", "not_a_probe_id"), 'hgu133plus2.db', link = 'bioGPS')
}
\seealso{
  \code{\link{biocann_orgdb}}
}

