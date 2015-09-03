\name{curated_expressiondata}
\alias{curated_expressiondata}
\alias{expression_curated}
\docType{data}
\title{
A toy dataset ofcurated RNA-seq to test within SMITE
}
\description{
A toy dataset of pre-cleaned gene expression data from RNA-seq.  The file is 
effect and p-value with gene names as rownames.
}
\usage{data("curated_expressiondata")}
\format{
  A data frame with 20819 observations on the following 2 variables.
  \describe{
    \item{rownames}{a character vector specifying gene}
    \item{column1}{an numeric vector specifying effect (log fold change)}
      \item{column2}{a numeric vector with a two sided p-value from DESeq 
      analysis}
  }
}
\value{
A dataframe with rownames as genes in Refseq format and columns for effects
and pvalues derived from negative binomial testing of DESeq normalized 
values from RNA-seq.
}
\details{
This gene expression dataset is a randomized version of the Toxoplasma dataset 
used to benchmark 
SMITE.  It no longer has NAs or  p-values=0.  Gene names were converted to gene 
symbols.
}
\source{
Manuscript in preparation.  Please see 
https://github.com/GreallyLab 
for more details.}
\examples{
data(curated_expressiondata)
}
\keyword{datasets}