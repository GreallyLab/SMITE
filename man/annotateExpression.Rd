\name{annotateExpression}
\alias{annotateExpression}
\alias{annotateExpression,PvalueAnnotation-method}
\title{
Adding expression data to a PvalueAnnotation
}
\description{
This function is used to create and load an ExpressionSet into a 
PvalueAnnotation.  Using specified effect and p-value column or named columns 
that the function will use to determine the effect and p-value columns, it loads 
the data it into the PvalueAnnotation. 
}
\usage{
annotateExpression(annotation, expdata, effect_col = NULL, pval_col = NULL)
}
\arguments{
  \item{annotation}{
An S4 object of class PvalueAnnotation
}
  \item{expdata}{
An object of class data.frame or matrix with row names corresponding to 
genes and atleast two columns with an effect and p-value for expression.

}
  \item{effect_col}{
A numeric specifying the column with an effect direction.  If not specified the 
function will grep for a single named column from: "effect", "odds",
"coeff" or "B"
}
  \item{pval_col}{
A numeric specifying the column with p-values.  If not specified the 
function will grep for a single named column from: "pval","p.val","p_val"
or "sig"
}
}
\details{
The function will load the entire given expression dataset as an ExpressionSet
in the expression slot,while the effect and p-value data will also be stored as 
an "AnnotatedDataFrame" in the phenoData slot of the ExpressionSet.
}
\value{
A PvalueAnnotation, an S4 object with the slot "expression" filled in.
}
\author{
N. Ari Wijetunga
}
\seealso{

annotateModification
makePvalueAnnotation
createPvalueObject
}
\examples{
##Note commented out. See vignette for more detailed usage information ##

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
\keyword{ expression }
\keyword{ SMITE }