\name{Reactome.Symbol.Igraph}
\alias{Reactome.Symbol.Igraph}
\alias{REACTOME}
\docType{data}
\title{
An Igraph network for REACTOME with nodes as gene symbols 
}
\description{
This is an Igraph network that was created using the REACTOME protein-protein 
interaction database.
}
\usage{data("Reactome.Symbol.Igraph")}
\format{
  An igraph object with 5770 nodes and 114288 edges
  \describe{
    \item{nodes}{gene names as gene symbols}
    \item{edges}{paired genes that interact}
  }
}
\details{
The provided igraph file was created using the igraph package and the 
interaction file provided from the reference. 
}
\value{
An Igraph network  based off of REACTOME interactions
}

\source{
REACTOME
}

\references{
http://www.reactome.org/pages/download-data/
}
\examples{
data(Reactome.Symbol.Igraph)
head(igraph::V(REACTOME))
}
\keyword{datasets}
