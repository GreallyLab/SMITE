\name{extractScores}
\alias{extractScores}
\alias{extractScores,PvalueAnnotation-method}
\title{
Extract scores for all genes
}
\description{
A function to obtain all gene scores
}
\usage{
extractScores(annotation)
}
\arguments{
  \item{annotation}{
An S4 object of class PvalueAnnotation for which scores have already
been calculated
}
}

\value{
A named vector containing all gene scores
}
\author{
N. Ari Wijetunga
}
\seealso{
scorePval
extractModules
highScores
}
\examples{
data(test_annotation_score_data)

out<-extractScores(test_annotation)

head(out)

}
