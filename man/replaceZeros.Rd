\name{replaceZeros}
\alias{replaceZeros}
\alias{replaceZeros,data.frame-method}
\title{
Replace Zeros in a Vector
}
\description{
This is a useful function to remove zeros and replace them with the minimum 
non-zero value of that vector
}
\usage{
replaceZeros(outdata, column)
}
\arguments{
  \item{outdata}{
 A dataframe
}
 \item{column}{
 A numeric indexing a particular column of the data frame
}
}
\details{
This function finds the smallest non-zero value in a particular column 
and then replaces all zeros with it within that column.  It will not do anything 
if there are no zeros.
}
\value{
The dataframe is returned.
}

\author{
N. Ari Wijetunga
}

\examples{

data (methylationdata)


methylation<-replaceZeros(methylation,5)

head(which(methylation[,5]==0))

}

\keyword{ zeros }
\keyword{ QC}
