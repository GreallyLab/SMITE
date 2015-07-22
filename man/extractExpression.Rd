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
extractExpression(annotation)
}

\arguments{
  \item{annotation}{
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
## Note: Commented out here.  See vignette for more detailed usage information #

## Load genome bed file ##
#data(hg19_genes_bed)

## Load curated expression bed file ##
#data(curated_expressiondata)

## Create a PvalueAnnotation with defaults for promoter size##
#test_annotation<-makePvalueAnnotation( data=hg19_genes, gene_name_col=5)

## Load Expression data into PvalueAnnotation ##
#test_annotation<-annotateExpression(test_annotation, expression_curated)

## Extract entire ExpressionSet with expression data ##
#slot(test_annotation,"expression")

## Extract expression data summary ##
#head(extractExpression(test_annotation))


}
\keyword{ SMITE }
\keyword{ expression }