\name{extractModSummary}
\alias{extractModSummary}
\alias{extractModules,PvalueAnnotation-method}
\title{
View summary of modification data within PvalueAnnotation
}
\description{
This function allows the user to view the metadata of the modification data
within the PvalueAnnotation object
}
\usage{
extractModSummary(pvalue_annotation)
}
\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation for which modification data has already
been loaded via annotateModification}
}
\value{
A data.frame pulled from the metadta slot of the modifications slot within a 
loaded PvalueAnnotation.
}
\author{
N. Ari Wijetunga
}
\seealso{
annotateModification
}
\examples{
##NOTE: Comment out in exmaple see vignette for more detailed usage ##

## Load genome bed file ##
data(hg19_genes_bed)

## Load curated DNA methylation bed file ##
#data(methylationdata)
#methylation<-methylation[-which(is.na(methylation[,5])),]
#methylation[,5]<-replace(methylation[,5],methylation[,5]==0, min(subset(methylation[,5], methylation[,5]!=0), na.rm=T))

## Create a PvalueAnnotation with defaults for promoter size##
test_annotation<-makePvalueAnnotation(data=hg19_genes, gene_name_col=5)

## Load DNA methylation into PvalueAnnotation ##
#test_annotation<-annotateModification(pvalue_annotation=test_annotation, 
#methylation, weight_by=c(promoter="distance", body="distance"), verbose=TRUE,
#mod_corr=FALSE, mod_type="methylation")

## Extract modification metadata ##
#extractModSummary(pvalue_annotation=test_annotation)
}

