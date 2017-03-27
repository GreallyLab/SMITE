\name{hg19_genes_bed}
\alias{hg19_genes_bed}
\alias{hg19_genes}
\docType{data}
\title{
A bed file annotating Refseq genes for the hg19 genome build
}
\description{
A gene anntation BED file containing columns for RefSeq name and Gene Symbol
}
\usage{data("hg19_genes_bed")}
\format{
  A data frame with 41633 observations on the following 6 variables.
  \describe{
    \item{column 1}{a character vector for chromosome}
    \item{column 2 }{an integer vector for start position}
    \item{column 3 }{an integer vector for end position}
    \item{column 4 }{an character vector for RefSeq gene names}
    \item{column 5 }{an character vector for Gene Symbol names}
    \item{column 6 }{an character vector for strand}
  }
}
\details{
A BED files taken from the table browser.
}
\value{
A dataframe in BED format (chromosome, start, end) with additional columns for
gene name as Refseq and gene symbol and strand.
}
\source{
Karolchik D, Hinrichs AS, Furey TS, Roskin KM, Sugnet CW, Haussler D, Kent WJ. 
The UCSC Table Browser data retrieval tool. Nucleic Acids Res. 2004 Jan 1;32
(Database issue):D493-6.
}
\references{
http://genome.ucsc.edu/
}
\examples{
data(hg19_genes_bed)
}
\keyword{datasets}
