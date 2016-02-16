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
annotateExpression(pvalue_annotation, expr_data, effect_col = NULL, pval_col = NULL)
}
\arguments{
    \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation
}
    \item{expr_data}{
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

data(curated_expressiondata)
data(test_annotation_score_data)
## Load Expression data into PvalueAnnotation ##
test_annotation <- annotateExpression(pvalue_annotation=test_annotation, 
expression_curated)

## Extract entire ExpressionSet with expression data ##
#slot(test_annotation,"expression")

## Extract expression data summary ##
#head(extractExpression(pvalue_annotation=test_annotation))

}
\keyword{ expression }
\keyword{ SMITE }