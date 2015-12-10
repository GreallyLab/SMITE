\name{plotCompareScores}
\alias{plotCompareScores}
\alias{plotCompareScores,PvalueAnnotation-method}
\title{
Compare two genomic features by score and display them in a hexbin plot
}
\description{
This function creates a hexbin of the log transformed p-value/score for any
two expression or modification-context pairing within a PvalueObject 
inside of a PvalueAnnotation
}
\usage{
plotCompareScores(pvalue_annotation, x_name, y_name, ...)
}
\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation for which a PvalueObject has already
been created
  }
  \item{x_name}{
A string to be grepped from the columns within the slot "pval_data" that is
within the PvalueAnnotation slot "score_data."  This column will be plotted on
the x-axis with a  direction specified from the corresponding effect column. 
}
  \item{y_name}{
A string to be grepped from the columns within the slot "pval_data" that is
within the PvalueAnnotation slot "score_data."  This column will be plotted on
the y-axis with a  direction specified from the corresponding effect column.
}
  \item{\dots}{
Other plotting parameters
}
}
\details{
This plotting function creates a hexbin plot of any two p-value vectors
stored within a p-value object.  It can be used to define relationships
between direction and significance in different genomic contexts after 
having combined p-values.
}
\value{
A hexbin plot
}
\author{
N. Ari Wijetunga
}
\seealso{
makePvalueObject
plotDensityPval
}
\examples{
data(test_annotation_score_data)

plotCompareScores(pvalue_annotation=test_annotation, x_name="expression", 
y_name="methylation_promoter")
}