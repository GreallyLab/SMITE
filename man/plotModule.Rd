\name{plotModule}
\alias{plotModule}
\alias{plotModule,PvalueAnnotation-method}
\title{
Plot a specific module after running Spinglass algorithm
}
\description{
This function is an adapted version of renderModule available through Epimods.
We have added optional functionality including plotting the actual raw data
onto the node edges, adding goseq annotation to the plot, legends and plotting 
modes.
}
\usage{
plotModule(pvalue_annotation, p_thresh = 0.05, which_network = 1, goseq = FALSE, 
layout = "fr", legend = TRUE, namestyle = "symbol", suppress_details = FALSE, 
meth_hi_col = "blue", meth_low_col = "yellow1", 
meth_mid_col = "gray90", exp_hi_col = "red1", exp_low_col = "chartreuse1", 
exp_mid_col = "gray90", label_scale = TRUE, compare_plot=FALSE, 
pdf_out=NULL)
}
\arguments{
  \item{pvalue_annotation}{
An S4 object of class PvalueAnnotation
}
  \item{p_thresh}{
A numeric specifying a threshold for plotting raw data on edges of
nodes.  DEFAULT is alpha=0.05.  Items above this threshold will be classified as
"mid" instead of "high" or "low"
"}
  \item{which_network}{
A numeric specifying which network to plot.  DEFAULTS to 1, and will
not plot another network until specified explicitly.
}
  \item{goseq}{
A logical indicating whether to plot goseq results for the module on the right 
hand side of the plot.
}
  \item{layout}{
A character string  as either "fr" (DEFAULT) for fruchterman.reingold or 
"circle" for a circular plot.
}
  \item{legend}{
A logical (TRUE(DEFAULT)/FALSE) specifying whether a legend should be drawn.
}
  \item{namestyle}{
A character string as either "symbol" (DEFAULT) for gene symbols or "refseq" for
RefSeq genes.  If modules were performed on RefSeq genes, then the function will
plot with gene symbols so that it is more useful.
}
  \item{suppress_details}{
A logical (TRUE(DEFAULT)/FALSE) indicating whether border raw data information 
should be plotted.
}
  \item{meth_hi_col}{
A color to be associated with signficant modification data with
positive effects
}
  \item{meth_low_col}{
A color to be associated with signficant modification data with
negative effects}
  \item{meth_mid_col}{
A color to be associated with non-signficant modification data 
}
  \item{exp_hi_col}{
A color to be associated with signficant expression data with
postive effects
}
  \item{exp_low_col}{
A color to be associated with signficant expression data with
negative effects}
  \item{exp_mid_col}{
A color to be associated with non-signficant expression data 
}
  \item{compare_plot}{
A logical (TRUE/FALSE(DEFAULT)) indicating whether two plots should be drawn 
side by side, one with raw data and one without
}
  \item{label_scale}{
A logical (TRUE/FALSE(DEFAULT)) indicating whether whether the node label should
be scaled with the node score
}
  \item{pdf_out}{
A string indicating a location to which the function should output a pdf.  If 
NULL (DEFAULT) then no pdf is made.
}

}

\value{
A plot of the module
}
\author{
N. Ari Wijetunga}
\seealso{
extractModules
}
\examples{
data(test_annotation_score_data)

#plotModule(pvalue_annotation=test_annotation, which_network=2)

#plotModule(pvalue_annotation=test_annotation, which_network=2,
#suppressDetails=TRUE)
}
