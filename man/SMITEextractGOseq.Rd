\name{SMITEextractGOseq}
\alias{SMITEextractGOseq}
\alias{SMITEsearchGOseq}
\alias{SMITEsearchGOseq,PvalueAnnotation-method}
\alias{SMITEextractGOseq,PvalueAnnotation-method}
\title{
 View the GOseq pathway analysis after having run Goseq, or search for a term.
}
\description{
Having defined at least one genomic module using SMITErunSpinglass or 
SMITErunBioNet, this function allows you to interrogate the enriched terms for a
specific module or combination of modules.
}
\usage{
SMITEextractGOseq(object, which.network = NULL)
SMITEsearchGOseq(object, searchstring, wholeword = FALSE)
}
\arguments{
  \item{object}{
An S4 object of class PvalueAnnotation, for which module-finding and GOseq 
analysis have already been performed
}
  \item{which.network}{
A numeric vector of a length of at least one, corresponding to a particular 
functional module specifically for the extract function.
}
 \item{searchstring}{
A character specfying a search string specifically for the search function.
}
  \item{wholeword}{
A logical (TRUE/FALSE) detemring whether the search string must be matched for
whole word specifically for the search function.
}
}
\details{
Goseq analysis is useful since it allows you to assess term/pathway enrichment
in a collection of genes, while adjusting for bias data.  Potential bias
can be from aspects like gene length or probe density that influence the
likelihood of finding a particular gene.  For more information please see
the goseq reference.
}
\value{
##Extract##
A list with eahc element matching the specified module.  Has columns 
identfying the term id, the over represented p-value, underrepresented 
p-value the total number in the category found in the module, the total
number in the category and a more descirptive pathway name. 


##Search##
A matrix with columns identfying the module name, module position/significance
, the specific enriched term, the rank of that term within all enriched terms 
and the total number of enriched terms.
}
\references{
Young MD, Wakefield MJ, Smyth GK and Oshlack A (2010). 
Gene ontology analysis for RNA-seq: accounting for selection bias.
Genome Biology, 11, pp. R14.
}
\author{
N. Ari  Wijetunga
}

\seealso{
SMITErunGOseq
SMITErunSpinglass
SMITErunBioNet
SMITEextractModules
SMITEplotModule
}
\examples{

## Commented out below.  See vignette for more details ##

##load sample data with only PvalueObject filled in##
#data(test_annotation_score_data)

## show goseq analysis for module 1 ##
#SMITEextractGOseq(test_annotation,1)

## show goseq analysis for module 1 and 2 ##
#SMITEextractGOseq(test_annotation,1:2)

## search for term ##

#SMITEsearchGOseq(test_annotation,"Cell cycle")
}

\keyword{ goseq }
\keyword{ SMITE }
