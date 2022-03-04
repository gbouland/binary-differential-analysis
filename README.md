# binary-differential-analysis
## Differential analysis of binarized single-cell RNA sequencing data captures biological variation 

Single-cell RNA sequencing data is characterized by a large number of zero counts, yet there is growing evidence that these zeros reflect biological variation rather than technical artifacts. We propose to use binarized expression profiles to identify the effects of biological variation in single-cell RNA sequencing data. Using 16 publicly available and simulated datasets, we show that a binarized representation of single-cell expression data accurately represents biological variation and reveals the relative abundance of transcripts more robustly than counts. BDA is available at https://github.com/gbouland/BDA

The scripts, functions and source data for the figures of the associated manuscript can be found here. Additionally, two vignettes have been made describing a DDA starting from a Seurat object en a raw count matrix

The paper has been published in NAR Genomics and Bioinformatics: https://doi.org/10.1093/nargab/lqab118

The results can interactively be explored here: http://insyprojects.ewi.tudelft.nl:5000/BinaryDifferentialAnalysis/



