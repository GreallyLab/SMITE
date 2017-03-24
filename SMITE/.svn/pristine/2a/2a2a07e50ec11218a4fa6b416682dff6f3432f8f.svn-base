\name{genes_for_conversiontest}
\alias{genes_for_conversiontest}
\alias{genes}
\docType{data}
\title{
A small set of RefSeq genes for converting
}
\description{
This toy dataset has 100 randomly selected RefSeq genes and can be used to test 
conversion functionality in SMITE
}
\usage{data("genes_for_conversiontest")}
\format{
  A data frame with 100 observations on the following 1 variables.
  \describe{
    \item{column1}{a character vector of RefSeqGene IDs}
  }
}
\value{
A dataframe with genes in Refseq format for conversion testing.
}

\examples{

data("genes_for_conversiontest")

genes[,1]<-convertGeneIds(gene_IDs=genes[,1], ID_type="refseq",
                          ID_convert_to="symbol")

}
\keyword{datasets}