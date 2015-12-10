\name{makePvalueAnnotation}
\alias{makePvalueAnnotation}
\alias{makePvalueAnnotation,ANY-method}
\title{
Initialize a PvalueAnnotation
}
\description{
This function initializes a PvalueAnnotation using a gene BED file and optional
BED files corresponding to interval datasets.  This is a necessary first step
in order to establish for each gene the gene promoter, body and associated
intervals.
}
\usage{
makePvalueAnnotation(data, other_data = NULL, other_tss_distance = 10000, 
    promoter_upstream_distance = 1000, promoter_downstream_distance = 1000, 
    strand_col = NULL, gene_name_col = NULL)
}
\arguments{
  \item{data}{
A required gene annotation BED file like that obtained from the UCSC Table 
Browser. At a minimum BED files must have the first three columns as 
(chromosome, start, end).  Additional required columns should correspond to the 
strand and gene name.The gene name needs to match the gene format desired
for the interaction network.  Duplicated gene names and associated gene
annotations are removed.

}
  \item{other_data}{
A list of BED files corresponding to each additional interval file to be 
associated with genes.  The function will use the other_d variable and the gene 
transcription start site (TSS) to find for each gene all intervals
within [tss-other_tss_distance, tss+other_tss_distance].
}
  \item{other_tss_distance}{
A vector specifying for each element of otherdata a distance from the gene
TSS to consider that interval as related to a gene.  If the length of other_d
does not match the length of the otherdata list, then the first value is used
for all datasets in the otherdata list. DEFAULTS to 10,000 base pairs.
}
  \item{promoter_upstream_distance}{
A numeric specifying how far upstream from the gene TSS is considered part of
the gene promoter. DEFAULTS to 1,000 base pairs.
}
  \item{promoter_downstream_distance}{
A numeric specifying how far downstream from the gene TSS is considered part of
the gene promoter.  Gene bodies subtract the promoter_downstream region.
DEFAULTS to 1,000 base pairs.
}
  \item{strand_col}{
A numeric specifying the column of the gene BED file (data) corresponding to 
the gene strand.  If this is not provided, the function will attempt to 
determine the strand.
}
  \item{gene_name_col}{
A numeric specifying the column of the gene BED file (data) corresponding to 
the gene name.
}
}
\details{
    The required only input file is the gene annotation BED file that should 
have (as all BED files) the chromosome, start and end  in columns 1,2 and 3, 
respectively.  Also, there should be a column for gene name and gene strand. 
The user needs to determine distance from the gene transcription start site that 
will define the gene promoter.  The gene body will then be calculated as the 
non-promoter overlapping sequence.  If optional BED files are given as
otherdata (e.g. transcription factor binding sites, histone modification peaks),
then the user will also decide a distance from the gene TSS to associate each
BED interval with a gene.  For a particular BED file, each genes may have more 
than one interval that falls within the desired range around a TSS.
    Unique gene names are required and the function  will automatically remove
duplicated genes.  We recommend deciding on an interaction network first and 
then loading a gene annotation BED file with the same names.  This will likely
necessitate allowing the function to pick one annotation of a gene, or pre-
processing using some criteria (e.g. longest transcript).
}
\value{
An S4 object of class PvalueAnnotation containing slots for an annotation
(GRangesList), an expression set, modifications (GRangesList), and a 
PvalueObject.
}
\author{
N. Ari Wijetunga
}

\seealso{
SMITE vignette
}
\examples{
## Note: Commented out below. See vignette for more detailed usage information##

## Load genome bed file ##
data(hg19_genes_bed)

## Create a PvalueAnnotation with defaults for promoter size##
test_annotation <- makePvalueAnnotation(data=hg19_genes, gene_name_col=5)

}
\keyword{ SMITE }
\keyword{ annotation }