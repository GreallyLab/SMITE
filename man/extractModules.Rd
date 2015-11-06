\name{extractModules}
\alias{extractModules}
\alias{extractModules,PvalueAnnotation-method}
\title{
View specific modules within a PvalueAnnotation
}
\description{
Having identfied modules within a Pvalue annotation, this function allows the user to display 1 or more of the module genes.
}
\usage{
extractModules(annotation, whichModule = NULL)
}
\arguments{
  \item{annotation}{
An S4 object of class PvalueAnnotation for which Spinglass or BioNet has already been run.
}
  \item{whichModule}{
A numeric vector specifying one or more module to display
}
}
\value{
A list with each element containing the requested modules
}
\author{
N. Ari Wijetunga
}
\seealso{
plotModule
runGOseq
extractGOseq
runSpinglass
runBioNet
}
\examples{
data(test_annotation_score_data)

extractModules(test_annotation, 1)

}