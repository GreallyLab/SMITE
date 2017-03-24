\name{extractExpression}
\alias{extractExpression}
\alias{extractExpression,PvalueAnnotation-method}
\title{
View the expression data stored in a PvalueAnnotation
}
\description{
This function allows the user to see the effect and p-value data that was loaded
into a PvalueAnnotation before performing downstream analysis.
}
\usage{
extractExpression(pvalue_annotation)
}

\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation for which expression data has already
been loaded via annotateExpression
}
}
\value{
A data.frame pulled from the phenoData of the expression slot within a load
PvalueAnnotation.  The phenoData specifically hold the effect and p-value
information.
}
\author{
N. Ari Wijetunga
}
\seealso{
annotateExpression
}
\examples{
data(test_annotation_score_data)
data(curated_expressiondata)
## Load Expression data into PvalueAnnotation ##
test_annotation<-annotateExpression(test_annotation, expression_curated)

## Extract entire ExpressionSet with expression data ##
#slot(test_annotation,"expression")

## Extract expression data summary ##
head(extractExpression(pvalue_annotation=test_annotation))


}
\keyword{ SMITE }
\keyword{ expression }