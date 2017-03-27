\name{plotDensityPval}
\alias{plotDensityPval}
\alias{plotDensityPval,PvalueAnnotation-method}
\title{
Plot the density of the combined scores stored in a PvalueObject
}
\description{
This function in called by the normalizePval function, but can also be called
by the user to visualize the relative densities of combined p-values (scores). 
}
\usage{
plotDensityPval(pvalue_annotation, ref = "expression_pvalue", ...)
}
\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation.
}
  \item{ref}{
A character specifying the name of the reference category. DEFAULT is 
"expression_pvalue"
}
\item{...}{
 Additional plotting arguments
 }
}
\value{
Plots a multidensity plot.
}
\author{
N. Ari Wijetunga
}
\examples{
## Load test annotation with only score data ##
data(test_annotation_score_data)

plotDensityPval(pvalue_annotation=test_annotation)

}
\keyword{ density }
\keyword{ SMITE }