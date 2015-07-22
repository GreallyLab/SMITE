# SMITE
This R package **SMITE**, Significance-based Modules Integrating the Transcriptome and Epigenome, allows reasearchers to score genes using a class of gene centric annotations associated with p-values from any prior statistical results. It builds on the previous framework of Epimods of using the Spin-glass algorithm implemented in the iGraph package with random sampling to find optimal weighted subnetworks ("modules") within a gene interaction network.  Identified modules are functionally annotated to aid in interpretation of results.

##Motivation
The critical goal in genomics research is the interpretation of multi-layered datasets that can provide more biological insight when considered together.  Instead of pairwise comparisons to determine overlapping genomic effects, integrated approaches provide a powerful, quick, and intuitive discovery platform.

##Usage
```{r}
library(SMITE)
browseVignettes("SMITE") #Please see our Tutorial for an in depth description of commands and usage of SMITE
```

##License
This package is free and open source software, licensed under GPL(>=2).
