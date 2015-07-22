\name{SMITEextractScores}
\alias{SMITEextractScores}
\alias{SMITEextractScores,PvalueAnnotation-method}
\title{
Extract scores for all genes
}
\description{
A function to obtain all gene scores
}
\usage{
SMITEextractScores(annotation)
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
SMITEscorePval
SMITEextractModules
SMITEhighScores
}
\examples{
data(test_annotation_score_data)

out<-SMITEextractScores(test_annotation)

head(out)

}
