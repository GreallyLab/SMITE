\name{runGOseq}
\alias{runGOseq}
\alias{runGOseq,PvalueAnnotation-method}
\title{
Run a GoSeq pathway analysis
}
\description{
This function allows pathway annotation of identified modules.
}
\usage{
runGOseq(pvalue_annotation, p_thresh = 0.05, coverage, type = "reactome")
}

\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation, for which module-finding has already 
been performed
}
  \item{p_thresh}{
A numeric specifying a threshold for signficance of FDR (q-values). DEFAULT is
alpha=0.05
}
  \item{coverage}{
A data.frame that is a bed file (chr  start stop) folowed by a gene name and
a numeric specifying the bias data (e.g. gene length or number of probes related
to gene)}
  \item{type}{
Either "kegg" to run KEGG analysis or "reactome" to run a REACTOME analysis
}
}
\details{
Goseq analysis is useful since it allows you to assess term/pathway enrichment
in a collection of genes, while adjusting for bias data.  Potential bias
can be from aspects like gene length or probe density that influence the
likelihood of finding a particular gene.  For more information please see
the goseq reference.

The function will compare all of the genes within a module to known pathways and
terms to find the terms that are most enriched within a module.  In this way,
this tool allows a reasearch to find a functional importance of a module.

We currently offer KEGG and REACTOME analysis, although additional pathway tools
may be added in the near future.

}
\value{
A PvalueAnnotation with goseq annotated modules.
}
\references{
Young MD, Wakefield MJ, Smyth GK and Oshlack A (2010). 
Gene ontology analysis for RNA-seq: accounting for selection bias. 
Genome Biology, 11, pp. R14.
}
\author{
Matthew D. Young myoung at wehi.edu.au
}
\note{
This is a wrapper function written by N. Ari Wijetunga for the package SMITE.
}

\seealso{
searchGOseq
extractGOseq
runSpinglass
runBioNet
extractModules
plotModule
}
\examples{
##load sample data with only PvalueObject filled in##
data(test_annotation_score_data)

## NOTE commented out in example.  Please see Vignette for better example ##
#test_annotation<-runGOseq(pvalue_annotation=test_annotation,
#coverage=read.table(
#system.file("extdata", "hg19_symbol_hpaii.sites.inbodyand2kbupstream.bed.gz", 
#package="SMITE"),stringsAsFactors=FALSE),type="kegg")


## search for a term ##
searchGOseq(test_annotation,"Cell cycle")

## show goseq analysis for module 1 ##
#extractGOseq(test_annotation,1)

}

\keyword{ goseq }
\keyword{ SMITE }