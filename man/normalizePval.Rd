\name{normalizePval}
\alias{normalizePval}
\alias{normalizePval,PvalueAnnotation-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function normalizes p-values (Scores) that are otherwise on different 
scales.
}
\description{
This function is a used to rescale compenent scores when distributions
have been altered.  There are two methods available.}
\usage{
normalizePval(pvalue_annotation, trans, ref = "expression_pvalue", method = "rescale")
}
\arguments{
  \item{pvalue_annotation}{
An S4 object of class p-value annotation
}
  \item{trans}{
A vector specifying a specific Box-cox power transformation to use. If not 
specified, the optimal transformation powers will be estimated. 
}
  \item{ref}{
A string that will be grepped from the names of the loaded expression data or
modification/context pairing.  All scores will be rescaled to match this 
reference's range.  The DEFAULT is expression.
}
  \item{method}{
A string of either

                    "Rescale"   DEFAULT
                    "rescale"
                                Performs a logit transform and rescales all
                                probabilites to the reference's logit 
                                transformed scale, then back-transforms

                    "Box-Cox"
                    "box-cox"
                    "boxcox" 
                    "Boxcox"
                                For each probility vector does a logit transform
                                and then iterates between 0.5 and 0.95 by 0.05
                                to determine the most similar transformation
                                to the logit transformed referece by a Wilcoxon-
                                test
}
}
\details{
Normalization may not be necessary but should improve some p-values from driving
the majority of downstream scores and modules solely because of the scale of 
their p-values.  All transformations are monotonic and are controlled for by
use of randomization prodecure downstream.
procedures downstream should }
\value{
An S4 object of class PvalueAnnotation with normalized p-values withing the
pval_data slot of the PvalueObject "score_data" slot

Plots densities of p-values before and after transform
}
\author{
N. Ari Wijetunga
}
\seealso{
makePvalueObject
scorePval
plotDensityPval
}
\examples{
data(test_annotation_score_data)

#test_annotation<-normalizePval(pvalue_annotation=test_annotation)

}