# differential-dropout-analysis
## Differential dropout analysis captures biological variation in single-cell RNA sequencing data

Single-cell RNA sequencing data is characterized by a large number of zero counts, yet there is growing evidence that these zeros reflect biological rather than technical artifacts. We propose differential dropout analysis (DDA), as an alternative to differential expression analysis (DEA), to identify the effects of biological variation in single-cell RNA sequencing data. Using 16 publicly available datasets, we show that dropout patterns are biological in nature and can assess the relative abundance of transcripts more robustly than counts.

The scripts, functions and source data for the figures of the associated manuscript can be found here. Additionally, two vignettes have been made describing a DDA starting from a Seurat object en a raw count matrix

A preprint of the manuscript can be found on bioRxiv: https://doi.org/10.1101/2021.02.01.429187 
The results can be interactively explored in a shiny app: https://gerard-bouland.shinyapps.io/differentialdropoutanalysis/
