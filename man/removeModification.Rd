\name{removeModification}
\alias{removeModification}
\alias{removeModification,PvalueAnnotation-method}
\title{
A function to "unload" a modification that has already been added.
}
\description{
After using the annotateModification function to load a modification into a
PValue annotation, you may wish to remove a modification or reannotate one,
which requires removing it first.
}
\usage{
removeModification(pvalue_annotation, mod_type = "methylation")
}
\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation
}
  \item{mod_type}{
A character string that identifies a type of modification within a 
PvalueAnnotation.  
}
}
\value{
An S4 object of class PvalueAnnotation
}
\author{
N. Ari Wijetunga
}
\seealso{
annotateModification
extractModification
extractModSummary
}
\examples{
##NOTE: Commented out in example. ##
## Please see vignette for more detailed usage information ##

## Load genome bed file ##

#data(hg19_genes_bed)

## Load curated DNA methylation bed file ##

data(methylationdata)
methylation <- methylation[-which(is.na(methylation[, 5])), ]
methylation[, 5] <- replace(methylation[, 5], methylation[, 5] == 0, 
    min(subset(methylation[, 5], methylation[, 5] != 0), na.rm=TRUE))

#meth1<-methylation

## make second curated test methylation bed file ##

#meth2<-methylation

## Create a PvalueAnnotation with defaults for promoter size##
#test_annotation<-makePvalueAnnotation(data=hg19_genes, gene_name_col=5)

## Load DNA methylation into PvalueAnnotation ##
#test_annotation<-annotateModification(annotation=test_annotation, 
#mod_data=meth1, weight_by=c(promoter="distance",body="distance"),verbose=TRUE,
#mod_corr=TRUE,mod_type="methylation")

## Extract GRanges with modification data ##
#extractModification(test_annotation)

## Load second dataset bed file ##
#test_annotation<-annotateModification(pvalue_annotation=test_annotation, 
#mod_data=meth2, weight_by=c(promoter="distance",body="distance"),
#verbose=TRUE, mod_corr=TRUE,mod_type="hydroxy")

## Extract GRanges with both modification dataset loaded ##
#head(extractModification(test_annotation,"hydroxy"))

## Unload DNA hydroxymethylation form PvalueAnnotation ##

#test_annotation<-removeModification(pvalue_annotation=test_annotation, 
                                    #mod_type="hydroxy")

## Extract GRanges to see only one modification dataset loaded ##
#head(extractModification(pvalue_annotation=test_annotation))

}

\keyword{ SMITE }
\keyword{ epigenetics }