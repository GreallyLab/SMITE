\name{methylationdata}
\alias{methylationdata}
\alias{methylation}
\docType{data}
\title{
A toy dataset of DNA methylation to test within SMITE
}
\description{
A toy dataset of raw DNA methylation from HELP-tagging.  The file is a 
BED file with columns added for effect and p-value.
}
\usage{data(methylationdata)}
\format{
  A data frame with 40000 observations on the following 5 variables.
  \describe{
    \item{column1}{a character vector specifying chromsome}
    \item{column2}{an integer vector specifying start}
    \item{column3}{an integer vector specifying end}
    \item{column4}{a numeric vector with an effect direction (here it is average
    difference between two groups)}
    \item{column5}{a numeric vector with a two sided t-test p-value}
  }
}
\details{
This is a small subset of a DNA methylation dataset is a randomized version of 
the Toxoplasma dataset used to benchmark.  We could not include the larger
version do to package size requirements but larger versions are available.
See Github source below. It still has NAs and p-values=0.
}
\source{
Manuscript in preparation. Please see https://github.com/GreallyLab/SMITE 
for more details.
}
\value{
A dataframe in BED format (chromosome, start, end) with additional columns for
and effect direction and p-value derived from T-tests of HELP-tagging DNA
methylation data.
}
\examples{
data(methylationdata)

any(is.na(methylation[, 4]))

any(methylation[, 4] == 0)
}
\keyword{datasets}