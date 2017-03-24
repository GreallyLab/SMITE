\name{runSpinglass}
\alias{runSpinglass}
\alias{runSpinglass,PvalueAnnotation-method}
\title{
Run Spinglass algorithm on a Scored PvalueAnnotation
}
\description{
This function is a  function to prepare the data for calling the Spinglass 
network algorithm.
}
\usage{
runSpinglass(pvalue_annotation, network, random_alpha = 0.05, gam = 0.5, 
node_alpha = 0.05, maxsize = 500, minsize = 8, num_iterations = 1000, simplify = TRUE)
}
\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation
}
  \item{network}{
An graph object of class graphNEL or igraph
}
  \item{random_alpha}{
A numeric specifying a threshold with with to determine module signficance 
after randomization
}
  \item{gam}{
A parameter used by the Spinglass algorithm
}
  \item{node_alpha}{
The proportion of nodes to be used as seeds for the community detection
}
  \item{maxsize}{
The maximum module size
}
  \item{minsize}{
The minimum module size
}
  \item{num_iterations}{
The number of randomizations that will be computed to determine whether the
module is significant by chance
}
 \item{simplify}{
A logical (TRUE(DEFAULT)/FALSE) that specifies whether network should be
simplified by removing self loops and repeated edges
}
}
\details{
In the provided Epimods reference, West et al outlined the advantages of using 
the spin-glass algorithm in the detection of modules.  Please consult the 
reference for more detailed information on the spin-glass algorithm  implemented
in the package igraph.

Like Epimods, this function employs the spin-glass algorithm implemented in 
igraph and uses random permutations to assess the "modularity," the number and 
strength of connected nodes, of a module.  However, SMITE scores are interpreted
as Chi-square distributed statistics whenever possible, rather than the
weighted-T-statistic in Epimods.

}
\value{
An S4 object of class PvalueAnnotation with modules loaded
}
\references{
James West, Stephan Beck, Xiangdong Wang & Andrew E. Teschendorff
An integrative network algorithm identifies age-associated differential 
methylation interactome hotspots targeting stem-cell differentiation pathway. 
Scientific Reports 3, Article number: 1630 (2013)

https://code.google.com/p/epimods/

}
\author{
N. Ari Wijetunga
}
\note{
This function was adapted from a function in the Epimods package that 
employs the spin-glass algorithm and uses random permutations to assess
the "modularity" of a module . The original function was created by West et al.
}

\seealso{
FEM
runBioNet
extractModules
plotModule
}
\examples{
data(test_annotation_score_data)

#load(system.file("data","Reactome.Symbol.Igraph.rda", package="SMITE"))
 
## NOTE: commented out for example.  See vignette for better explanation ##
#test_annotation <- runSpinglass(pvalue_annotation=test_annotation, 
#network=REACTOME, maxsize=50, num_iterations=10)

plotModule(test_annotation, which_network=6, layout="fr")

}
