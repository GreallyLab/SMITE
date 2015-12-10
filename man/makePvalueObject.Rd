\name{makePvalueObject}
\alias{makePvalueObject,PvalueAnnotation-method}
\alias{makePvalueObject}
\title{
Function to make a PvalueObject within a PvalueAnnotation}
\description{
Having annotated modifications and expression data this function will assemble
a PvalueObject within the slot "score_data" of a PvalueAnnotation.  This
is a necessary step before being able to run downstrem functions.
}
\usage{
makePvalueObject(pvalue_annotation, effect_directions = NULL)
}
\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation
}
  \item{effect_directions}{
A character vector with optional names specifying
                "increase"  Modification is expected to increase as expression
                            increase
                "decrease"  Modification is expected to descrease as 
                            expression decreases
                "bidirectional"
                            No direction is assumed between modification and
                            direction
}
}
\details{
The specified relationship between the modification and expression will be 
stored and then used when scoring.
}
\value{
An S4 object of class PvalueAnnotation with a slot for score_data filled it
}
\author{
N.Ari Wijetunga
}

\seealso{
makePvalueAnnotation
}
\examples{
#NOTE: Commented out in example, please see vignette for more details##
options(stringsAsFactors=FALSE)

data(methylationdata)
methylation <- methylation[-which(is.na(methylation[, 5])), ]
#methylation[, 5] <- replace(methylation[, 5],methylation[, 5] == 0, 
#min(subset(methylation[, 5], methylation[, 5] != 0), na.rm=TRUE))

#data(curated_expressiondata)
#data(hg19_genes_bed)
#data(histone_h3k4me1)

#test_annotation<-makePvalueAnnotation(data=hg19_genes, 
#other_data=list(h3k4me1=h3k4me1), gene_name_col=5,other_tss_distance=5000)

#fill in expression data
#test_annotation<-annotateExpression(test_annotation, expression_curated)


#fill in methylation data
#this step takes ~10 minutes
#test_annotation<-annotateModification(test_annotation, methylation, 
#weight_by=c(promoter="distance",body="distance",h3k4me1="distance"),
#verbose=TRUE, mod_corr=FALSE)

#create a pvalue object that will count the effect of the h3k4me1 as
#bidirectional
#test_annotation<-makePvalueObject(pvalue_annotation=test_annotation, 
#effect_directions=c(methylation_promoter="decrease",
#methylation_body="decrease", methylation_h3k4me1="bidirectional"))
}
\keyword{ SMITE }