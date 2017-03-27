\name{highScores}
\alias{highScores}
\alias{highScores,PvalueAnnotation-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a vector of the highest scoring genes
}
\description{
This function can be used to extract a subset of the highest scoring genes for
other downstream analysis. 
}
\usage{
highScores(pvalue_annotation, alpha = 0.05)
}
\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation for which scoring has already been
performed
}
  \item{alpha}{
A numeric specifying a threshold at signficant genes can be determined. DEFAULT
is alpha=0.05.
}
}
\details{
This function randomly samples the scores with replacement 100 times and for
within each random sample for each score it determines the proportion of scores
at or greater than the score. The average of these proportions over the 100 
samples will be the new p-value/scores.  All scores falling below the threshold
will be returned.

}
\value{
A named vector of scores.
}
\author{
N. Ari Wijetunga
}
\seealso{
scorePval
plotCompareScores
runSpinglass
runBioNet
}
\examples{
data(test_annotation_score_data)


## Note: commented out for example.  See vignette for better example ##

#out <- highScores(pvalue_annotation=test_annotation, alpha=0.01)

}
\keyword{ nonparametric }
\keyword{ genes }
