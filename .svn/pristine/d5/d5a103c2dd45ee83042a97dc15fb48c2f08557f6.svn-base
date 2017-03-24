\name{SMITE-package}
\alias{SMITE-package}
\alias{SMITE}
\docType{package}
\title{
Significance-based Modules Integrating the Transcriptome and Epigenome
}
\description{
SMITE provides a method of scoring and visualizing multi-level epigenomic data 
in order to prioritize genes within a genome-wide experiment.  These scores can 
then be used to identify subnetworks within an interaction network
called modules.  Each module represents a collection of highly interacting genes
that are implicated by the experiment.
}
\details{
\tabular{ll}{
Package: \tab SMITE\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2015-07-06\cr
License: \tab GPL (>=2)\cr
}
}
\author{
Neil Ari Wijetunga,
Andrew Damon Johnston

Maintainer: Neil.Wijetunga@med.einstein.yu.edu, Andrew.Johnston@med.einstein.yu.edu
}

\keyword{ package }
\seealso{
FEM
BioNet
}
\examples{
## NOTE: commented out for example.  See vignette for better explanation ##

options(stringsAsFactors=FALSE)

data(methylationdata)
methylation <- methylation[-which(is.na(methylation[, 5])), ]
methylation[,5] <- replace(methylation[,5],methylation[,5] == 0, 
    min(subset(methylation[,5], methylation[,5]!=0), na.rm=TRUE))

data(curated_expressiondata)

data(hg19_genes_bed)
data(histone_h3k4me1)

#test_annotation<-makePvalueAnnotation( data=hg19_genes, 
#other_data=list(h3k4me1=h3k4me1), gene_name_col=5, other_tss_distance=5000)

##fill in expression data
#test_annotation<-annotateExpression(test_annotation, expression_curated)


##fill in methylation data

#test_annotation<-annotateModification(test_annotation, methylation, 
#weight_by=c(promoter="distance", body="distance", h3k4me1="distance"),
#verbose=TRUE, mod_corr=TRUE)

##create a pvalue object that will count the effect of the h3k4me1 as
##bidirectional

#test_annotation<-makePvalueObject(test_annotation, 
#effect_directions=c(methylation_promoter="decrease",
#methylation_body="decrease",
#methylation_h3k4me1="bidirectional"))

##normalize the pvalues compared to colExp

#test_annotation<-normalizePval(test_annotation,ref="expression_pvalue",
#method="rescale")
    
##score with all four features contributing

#test_annotation<-SMITEscorePval(test_annotation,
#weights=c(methylation_promoter=.3,methylation_body=.1,expression=.3,
#methylation_h3k4me1=.3))

##load REACTOME 
#load(system.file("data","Reactome.Symbol.Igraph.rda", package="SMITE"))
 
##run Spinglass using REACTOME network

#test_annotation<-runSpinglass(test_annotation, REACTOME, maxsize=50, 
#num_iterations=10)

##run goseq on individual modules to determine bias 
#test_annotation <- runGOseq(test_annotation,
#coverage=read.table(system.file("extdata", 
#"hg19_symbol_hpaii.sites.inbodyand2kbupstream.bed.gz", package="SMITE")),
#type="kegg")

##search go seq output for keywords
#searchGOseq(test_annotation, "Cell")

##Draw a network
#plotModule(test_annotation, which_network=6, layout="fr")

##sample final file ##
data(test_annotation_score_data)


}
