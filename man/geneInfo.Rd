\name{geneInfo}
\alias{geneInfo}
\title{Simple Feature Annotation}
\usage{
  geneInfo(x, annotation = "human",
    extras = c("biogps", "ncbi", "pathway", "kegg"))
}
\arguments{
  \item{x}{feature identifiers to annotate}

  \item{annotation}{Annotation package to use.
It can be the name of an annotation package or of an organism supported by
\code{\link{biocann_orgdb}}.}

  \item{extras}{Indicates the type of information/resource to add for each feature.
Out-links to online resource can be added the prefix \code{'~'} to selected resources.

Note that the resources \code{'NCBI'} and \code{'bioGPS'} are \emph{de facto} online resources,
and make ENTREZID and Symbol as links to the respective gene page.}
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
geneInfo(1:5, extras = 'bioGPS')

# probe ids
ids <- c("1007_s_at", "1053_at", "117_at", "121_at", "not_a_probe_id")
geneInfo(ids, 'hgu133plus2.db')
geneInfo(ids, 'hgu133plus2.db', extras = 'bioGPS')
geneInfo(ids, 'hgu133plus2.db', extras = 'pathwayID')
}
\seealso{
  \code{\link{biocann_orgdb}}
}

