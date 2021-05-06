library(Seurat)
library(magrittr)

meta <- rio::import("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/Mouse_aging/GSM4505405_tabula-muris-senis-facs-official-raw-obj-metadata.csv.gz")
all <- ReadH5AD("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/Mouse_aging/GSM4505405_tabula-muris-senis-facs-official-raw-obj.h5ad")
counts <- all@assays$RNA@counts
samplesheet <- meta[,c("index","tissue","age","mouse.id")]
samplesheet <- samplesheet[samplesheet$age %in% c("3m","24m"),]
tissues <- table(samplesheet$tissue,samplesheet$age) %>% as.matrix()
tissues <- tissues[tissues[,1]!= 0 & tissues[,2]!= 0,] %>% rownames()
samplesheet <- samplesheet[samplesheet$tissue %in% tissues,]
colnames(samplesheet) <- c("IDs","celltype","status","patient")
counts <- counts[,samplesheet$IDs]
Aging_mouse_FACS <- list("counts" = counts, "samplesheet" = samplesheet)
saveRDS(object = Aging_mouse_FACS, file = "C:/Users/gabouland/Documents/004 PhD/008scBinaryAnalysis/000 Preprocess/Aging_mouse_FACS.rds")


