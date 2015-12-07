\name{scorePval}
\alias{scorePval}
\alias{scorePval,PvalueAnnotation-method}
\title{
Making a single combined score for each gene
}
\description{
This function uses an a priori weighting scheme to combine scores for a given 
gene. 
}
\usage{
scorePval(pvalue_annotation, weights)
}
\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation, for which makePvalueObject has
already been run.
}
  \item{weights}{
A numeric vector of the relative importance of expression, modifications,
and genomic contexts toward the final score.  Names should be provided that 
match the "modfication_genomicfeature" format, except for expression.  While
the scores do not have to add up to 1, it is good practice to impose this 
restriction in order to track the relative importance.
}
}
\details{
Because each weighting scheme generates scores from a distribution that
will change depending on the analysis inputs, the function will randomly sample 
the final scores and compare each derived score to this simulated distribution.

If no names are given, then the function will assume the weights are in the
order that it finds a particular "modification_genomicfeature" and it will 
print the weighting scheme so that you can verify it is correct.  The total 
number of weights must match the total number of modifications*genomicfeatures+1
for expression.

After calculating a combined score (using a Stouffer's weighted statistic), 
a  new p-value is derived using a non-parametric sampling approach.
}
\value{
An S4 object of class PvalueAnnotation.
}
\author{
N. Ari Wijetunga
}
\examples{
options(stringsAsFactors=FALSE)

data(test_annotation_score_data)

## NOTE: commented out for example.  See vignette for better explanation ##
#test_annotation<-scorePval(pvalue_annotation=test_annotation,
#weights=c(methylation_promoter=.3,methylation_body=.1,expression=.3,
#methylation_h3k4me1=.3))

}

\keyword{ MonteCarlo }
\keyword{ SMITE }