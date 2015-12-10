\name{extractModification}
\alias{extractModification}
\alias{extractModification,PvalueAnnotation-method}
\alias{extractModSummary}
\alias{extractModSummary,PvalueAnnotation-method}
\title{
Extract some or all loaded modifications or a the summary of combined effects
}
\description{
After having loaded modifications into a PvalueAnnotation, these functions can 
be used to display the GRanges with the modification of interest, or the data
frame containg a summary of the combined effects.
}
\usage{
extractModification(pvalue_annotation, mod_type = "methylation")
extractModSummary(pvalue_annotation)
}

\arguments{
  \item{pvalue_annotation}{
An s4 object of calss PvalueAnnotation
}
  \item{mod_type}{
A string or character vector that must match one or more of the loaded 
modifications. If NULL (DEFAULT) then it will show all modifications.
}
}
\value{
A GRanges object containing the modification(s) of interest or a data frame with
a summary of the combined p-values and effects
}
\author{
N. Ari Wijetunga
}
\seealso{
extractExpression
annotateModificaiton
removeModification
}
\examples{
##NOTE: Comment out in exmaple see vignette for more detailed usage ##

## Load genome bed file ##
data(hg19_genes_bed)

## Load curated DNA methylation bed file ##
#data(methylationdata)
#methylation <- methylation[-which(is.na(methylation[,5])),]
#methylation[, 5] <- replace(methylation[,5],methylation[,5] == 0, 
#min(subset(methylation[, 5], methylation[, 5] !=0 ), na.rm=TRUE))

## Create a PvalueAnnotation with defaults for promoter size##
test_annotation<-makePvalueAnnotation(data=hg19_genes, gene_name_col=5)

## Load DNA methylation into PvalueAnnotation ##
#test_annotation <- annotateModification(pvalue_annotation=test_annotation, 
#methylation, weight_by=c(promoter="distance", body="distance"), verbose=TRUE,
#mod_corr=FALSE, mod_type="methylation")

## Extract GRanges with modification data ##
#extractModification(pvalue_annotation=test_annotation)

}
\keyword{ SMITE }
\keyword{ epigenetics }
\keyword{ networks }
