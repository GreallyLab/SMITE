\name{convertGeneIds}
\alias{convertGeneIds,character,character,character-method}
\alias{convertGeneIds}
\title{
Convert between gene ids
}
\description{
A convenient function used to convert between gene ids from different
gene annotations.
}
\usage{
convertGeneIds(a, a_type, b_type, delim = NULL, verbose = FALSE)
}
\arguments{
  \item{a}{
A vector of gene names.
}
  \item{a_type}{
A character specifying the type of given annotation.  Currently one of "refseq",
"ensembleprot","uniprot" or "ensemble". In the case of b_type="entrez", "symbol"
}
  \item{b_type}{
A character specifying the type of desired annotation.  Currently one of 
"symbol" or in they case of a_type="symbol", "entrez"
}
  \item{delim}{
An optional character that will be removed from the beginning of each gene 
name.  It can be a long string.
}
  \item{verbose}{
TRUE/FALSE Should the function be verbose?  DEFAULTS to FALSE.
}
}
\details{
This is a very usefule function to efficently convert between gene ids.  It
currently relies on enumeration of each possible conversion, which has limited
it's use to mainly converting to gene symbol.

}
\value{
A character vector formatted to b_type
}
\author{
N. Ari Wijetunga < Neil.Wijetunga@med.einstein.yu.edu >
}
\note{
The function has enumerated combinations using AnnotationDBI.  We can provide 
additional functionality if needed.
}

\seealso{
AnnotationDBI, Biomart
}
\examples{

data(genes_for_conversiontest)

genes[,1]<-convertGeneIds(a=genes[,1], a_type="refseq",
                                     b_type="symbol")
                                                 
}
\keyword{SMITE }
\keyword{ convert }
\keyword{ ids }