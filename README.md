# differential-dropout-analysis
## Differential dropout analysis captures biological variation in single-cell RNA sequencing data

Single-cell RNA sequencing data is characterized by a large number of zero counts, yet there is growing evidence that these zeros reflect biological variation rather than technical artifacts. We propose differential dropout analysis (DDA) to identify the effects of biological variation in single-cell RNA sequencing data. Using 16 publicly available and simulated datasets, we show that DDA accurately detects biological variation and can assess the relative abundance of transcripts more robustly than methods relying on counts. DDA is available at https://github.com/gbouland/DDA

The scripts, functions and source data for the figures of the associated manuscript can be found here. Additionally, two vignettes have been made describing a DDA starting from a Seurat object en a raw count matrix

A preprint of the manuscript can be found on bioRxiv: https://doi.org/10.1101/2021.02.01.429187 

The results can be interactively explored in a shiny app: http://insyprojects.ewi.tudelft.nl:5000/DifferentialDropoutAnalysis/


