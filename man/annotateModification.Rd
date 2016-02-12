\name{annotateModification}
\alias{annotateModification}
\alias{annotateModification,PvalueAnnotation-method}
\title{
Adding modification data to a PvalueAnnotation
}
\description{
This function is the main "workhorse" function for SMITE because given a 
specific epigenetic modification (e.g. DNA methylation) it will 1) assess
an internal correlation structure and 2) aggregate the modification over all 
intervals associated with a gene in the "makePvalueAnnotation" function.
}
\usage{
annotateModification(pvalue_annotation, mod_data, weight_by = NULL, 
weight_by_method = "Stouffer", mod_included = NULL, mod_corr = TRUE, 
mod_type = "methylation", verbose = FALSE)
}
\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation.
}
  \item{mod_data}{
A dataframe or matrix derived from a bed file with thhe the first three columns
as (chromosome, start, end), column 4 is the effect, and column 5 is the
p-value.
}
  \item{weight_by}{
A vector with named elements specifying how modifications should be weighted 
within an interval.  Must be one of:
    
    "distance"  Use the distance from the gene TSS to weight the p-values
                and the combined effect such that events closer to the
                TSS are weighted more.  Log distances are used.
    
    "pval" 
    "p.value"
    "pvalue"
    "p_val"     (DEFAULT) Do not weight p-values but weight the combined 
                effect such by the signficance of that effect.
     
     ELSE       Do not weight p-values or the combined effect.
    
}
  \item{weight_by_method}{
A character specifying which method should be used to combine 
p-values.  Must be one of:

    "Stouffer"  (DEFAULT)
                Stouffer's method for combing pvalues involves first taking the 
                inverse standard normal
                CDF transformation of a vector of p-values followed by a 
                weighted sum
                creating a new Z score with a standard normal distribution
    "Fisher"
    "fisher"
    "chisq"
    "chi"
                Fisher's method involves summing the -2ln(p) for each of k p's
                which follows an approximate
                chi square distribution with 2k degrees of freedom
    "Sidak"
    "sidak"
    "minimum"
                Sidak's adjustment is essentially the minimum p-value, with an 
                added transformation to account for
                multiple comparisons.
    "binomial"
                The binomial probability assesses the probability of observing 
                the observed number of p-value below a threshold (alpha=0.05) 
                given the total number of p values and the probability of a 
                false positive.
}
  \item{mod_included}{
A vector of named elements specifying for which intervals in the annotation the
function should find combined scores (e.g. promoters).  If not specified the
assumption is that all type of intervals associated with a gene should be 
included.
}
  \item{mod_corr}{
A logical (TRUE/FALSE) specifying whether a correlation matrix should be
estimated.  The DEFAULT is TRUE.
}
  \item{mod_type}{
A character naming the modification that is being loaded.  The DEFAULT is
"methylation" and any modType string can be used, but will be referred to 
in downstream analysis.  A unique name must be used for each modification that 
is loaded.  When picking a variable modType should also avoid using "_" as it is
used to split column names containing modType.
}
  \item{verbose}{
A logical specifying if the user wants updates about the progress of the
function.
}
}
\details{
This function is the main "workhorse" function for SMITE because given a 
specific epigenetic modification (e.g. DNA methylation) it will 1) assess
an internal correlation structure and 2) aggregate the modification over all 
intervals associated with a gene in the "makePvalueAnnotation" function.

}
\value{
A an S4 object of class PvalueAnnotation with the slot modification 
(a GrangesList) filled in for each additional modification.
}
\references{
Fisher R. Statistical methods for research workers. Oliver and Boyd; 
Edinburgh: 1932.

Stouffer S, DeVinney L, Suchmen E. The American soldier: 
Adjustment during army life. Vol. 1. Princeton University Press; 
Princeton, US: 1949.

Sidak, Z. (1967). Rectangular confidence regions for the means of
multivariate normal distributions, Journal of the American Statistical 
Association 62, 626 633.
}
\author{
N. Ari Wijetunga
}
\seealso{
removeModification
annotateExpression
makePvalueAnnotation
createPvalueObject
}
\examples{
options(stringsAsFactors=FALSE)

## Commented out below See vignette for more detailed usage information ##
## Load genome bed file ##
#data(hg19_genes_bed)

## Create a PvalueAnnotation with defaults for promoter size##
#test_annotation<-makePvalueAnnotation(data=hg19_genes, gene_name_col=5)

## Load DNA methylation bed file ##
#data(methylationdata)
#methylation<-methylation[-which(is.na(methylation[,5])),]
#methylation[,5]<-replace(methylation[,5],methylation[, 5] == 0, 
#min(subset(methylation[,5], methylation[,5]!=0), na.rm=TRUE))

## Load DNA methylation into PvalueAnnotation modCorr=F for example##
## NOTE: Commented out below.  See vignette for better example ##
#test_annotation <- annotateModification(pvalue_annotation=test_annotation, 
#mod_data=methylation, weight_by=c(promoter="distance", body="distance"),
#verbose=FALSE, mod_corr=FALSE, mod_type="methylation")

}
\keyword{ SMITE }
\keyword{ epigenetic }