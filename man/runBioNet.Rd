\name{runBioNet}
\alias{runBioNet}
\alias{runBioNet,PvalueAnnotation-method}
\title{
Perform BioNet Analysis on a PvalueAnnotation
}
\description{
With BioNet, a researcher can find a single interconnected gene module using
the highest scoring genes generated in a PvalueAnnotation.  This function will
load the module into the PvalueAnnotation for visualization and downstream 
analysis.}
\usage{
runBioNet(pvalue_annotation, network, alpha = 0.05)
}
\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation that has already had scores generated.
}
  \item{network}{
An interaction network of class graphNEL or igraph.
}
  \item{alpha}{
A numeric specfiying a cutoff for high scoring genes to be return with the
highScores function.
}
}
\details{
The input of p-values to BioNet discussed in the BioNet vignette involves
first modeling p-values as a Beta-uniform mixture model to obtain the actual
corresponding probability function values.  Since our scoring method produces
p-values/scores that are uniform in distribution, we input them directly into
the BioNet algorithm.  For more details on BioNet see the reference or
runFastHeinz in the BioNet package.
}
\value{
A PvalueAnnotation with a loaded module.
}
\references{
Beisser et al.
BioNet: an R-Package for the functional analysis of biological networks.
Bioinformatics. 2010 Apr 15;26(8):1129-30. doi: 10.1093/bioinformatics/btq089. 
Epub 2010 Feb 25.

}
\author{
N. Ari Wijetunga
}
\note{
This is a wrapper function to run BioNet.  The actual BioNet 
code was created by Beisser et al.
}
\seealso{
plotModule
extractModule
runGOseq
}
\examples{
## load test data ##
data(test_annotation_score_data)

 
## NOTE: commented out for example.  See vignette for better explanation ##

#load reactome network with gene symbols ##
#load(system.file("data","Reactome.Symbol.Igraph.rda", package="SMITE"))

## run BioNet ##
#test_annotation<-runBioNet(pvalue_annotation=test_annotation,
#network = REACTOME)

## view module ##
#extractModules(pvalue_annotation=test_annotation, 1)

## plot module ##
#plotModule(pvalue_annotation=test_annotation, which.network=1)

}
\keyword{ SMITE }
\keyword{ network }