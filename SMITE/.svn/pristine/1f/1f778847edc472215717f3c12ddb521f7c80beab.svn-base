\name{test_annotation_score_data}
\alias{test_annotation_score_data}
\alias{test_annotation}
\docType{data}
\title{
A toy PvalueAnnotation
}
\description{
 This Pvalue annotation has only scoring data filled in to use in late
pipeline "SMITE" functions.  It can be used to skip the loading data phase
of analysis and test latter functionality.
}
\usage{data("test_annotation_score_data")}
\format{
  A PvalueAnnotation with the following slots
  \describe{
    \item{score_data}{a PvalueObject with slots corresponding to pval_data,
    effect_data, genes, signs_index, scores, trans, scoring_vector, and 
    module_otuput}
  }
}
\details{
This is a PvalueAnnotation which has had all of the pre-scoring data 
removed so that it is only usefule for using functions beginning with SMITE
and SMITE plotting functions.
}
\value{
A PvalueAnnotation with the score_data slot containing toy scores
}

\examples{
data(test_annotation_score_data)
plotDensityPval(test_annotation)
head(extractScores(test_annotation))

}
\keyword{datasets}