\name{stoufferTest}
\alias{stoufferTest}
\alias{stoufferTest,vector-method}
\title{
Stouffer's Test
}
\description{
This function performs a weighted Stouffer's mehtod of combining p-values.
}
\usage{
stoufferTest(pvalues, weights)
}
\arguments{
  \item{pvalues}{
A vector of p-values.
}
  \item{weights}{
Optional weights used when combining probabilites.  If no weights are given then the p-values are equally weighted.
}
}
\details{
For each p-value the inverse standard normal CDF is applied and  Z scores 
are derived. Z-scores are then summed and a new Z score is transformed back to 
a p-value.
}
\value{
A numeric p-value that represents the standard normal CDF of the combined 
Z statistic.
}
\references{

https://en.wikipedia.org/wiki/Fisher's_method

Stouffer S, DeVinney L, Suchmen E. The American soldier: Adjustment during army life. Vol. 1. Princeton University Press; Princeton, US: 1949.

}
\note{
This function was adapted from the function written on the Fisher's Method
wikipedia page.
}
\examples{
## Generate test weights ##
weights<-runif(10, 1,100)
weights<-sort(weights)

## Generate test p-values##
pvals<-runif(10,0,1)

## run stoufferTest ##
stoufferTest(pvalues = pvals, weights=1/weights)

}
\keyword{ pvalues }
\keyword{ significance }
