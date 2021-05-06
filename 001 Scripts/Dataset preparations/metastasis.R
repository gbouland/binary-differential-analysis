library(magrittr)
library(Seurat)

setwd("/tudelft.net/staff-umbrella/pQTL/000scBinary/000 Data/004 Metastasis")
samplesheet <- rio::import("GSE131907_Lung_Cancer_cell_annotation.txt")
samplesheet <- samplesheet[samplesheet$Sample_Origin %in% c("mLN","nLN","nLung","tLung"),]
samplesheet$tissue <- ifelse(grepl("Lung",samplesheet$Sample_Origin),"Lung","LN")
samplesheet$status <- ifelse(grepl("n",samplesheet$Sample_Origin),"normal","cancer")
samplesheet$status <- ifelse(grepl("t",samplesheet$Sample_Origin),"cancer",samplesheet$status)
table(samplesheet$tissue,samplesheet$Cell_type)
rawCounts <- rio::import("GSE131907_Lung_Cancer_raw_UMI_matrix.txt")
rownames(rawCounts) <- rawCounts$Index
rawCounts$Index <- NULL
rawCounts <- rawCounts[,samplesheet$Index]

## Lung ##
samplesheet_lung <- samplesheet[samplesheet$tissue == "Lung",]
rawCounts_lung <- rawCounts[,samplesheet_lung$Index]
samplesheet_lung <- samplesheet_lung[,c("Index","Cell_type","status","Sample")]
colnames(samplesheet_lung) <- c("IDs","celltype","status","patient")
Metastatis_Lung_set <- list("counts" = rawCounts_lung, "samplesheet" = samplesheet_lung)
saveRDS(object = Metastatis_Lung_set, file = "Metastatis_Lung.rds")

## LN ##

samplesheet_LN <- samplesheet[samplesheet$tissue == "LN",]
samplesheet_LN <- samplesheet_LN[samplesheet_LN$Cell_type %in% c("B lymphocytes","Myeloid cells","T lymphocytes"),]
rawCounts_LN <- rawCounts[,samplesheet_LN$Index]
table(samplesheet_LN$Cell_type,samplesheet_LN$status)


samplesheet_LN <- samplesheet_LN[,c("Index","Cell_type","status","Sample")]
colnames(samplesheet_LN) <- c("IDs","celltype","status","patient")
Metastatis_LN_set <- list("counts" = rawCounts_LN, "samplesheet" = samplesheet_LN)
saveRDS(object = Metastatis_LN_set, file = "Metastatis_LN.rds")
