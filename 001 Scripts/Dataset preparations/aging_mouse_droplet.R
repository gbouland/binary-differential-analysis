library(Seurat)
library(magrittr)

meta <- rio::import("./GSM4505404_tabula-muris-senis-droplet-official-raw-obj-metadata.csv")
all <- ReadH5AD("./GSM4505404_tabula-muris-senis-droplet-official-raw-obj.h5ad")
counts <- all@assays$RNA@counts

samplesheet <- meta[,c("index","tissue","age","mouse.id")]
samplesheet <- samplesheet[samplesheet$age %in% c("3m","24m"),]

tissues <- table(samplesheet$tissue,samplesheet$age) %>% as.matrix()
tissues <- tissues[tissues[,1]!= 0 & tissues[,2]!= 0,] %>% rownames()

colnames(samplesheet) <- c("IDs","celltype","status","patient")
counts <- counts[,samplesheet$IDs]

Aging_mouse_droplet <- list("counts" = counts, "samplesheet" = samplesheet)
saveRDS(object = Aging_mouse_droplet, file = "./Aging_mouse_droplet.rds")

