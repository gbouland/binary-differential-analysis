library(Seurat)
library(Matrix)
atlas <- readRDS("./cancerAtlas.RDS")
samplesheet <- atlas@meta.data
counts <- atlas@assays$RNA
rm(atlas)

##cancerAtlas_naive_vs_transitional_tcells##
contrast <- c("Naive-memory CD4 T cells", "Transitional memory CD4 T cells")
samplesheet <- samplesheet[samplesheet$cell_type %in% contrast,]
counts <- counts[,rownames(samplesheet)]
out <- rowSums(counts)
out <- names(out[out==0])
counts <- counts[!rownames(counts) %in% out,]
samplesheet <- data.frame("IDs" = rownames(samplesheet),
                          "celltype" = samplesheet$source,
                          "status" = samplesheet$cell_type,
                          "patient" = samplesheet$patient)

cancerAtlas_naive_vs_transitional_tcells <- list("counts" = counts, "samplesheet" = samplesheet)
saveRDS(object = cancerAtlas_naive_vs_transitional_tcells, file = "./cancerAtlas_naive_vs_transitional_tcells.rds")

##Regulatory T cells vs T helper cells
table(samplesheet$source,samplesheet$cell_type)
contrast <- c("Regulatory T cells", "T helper cells")
samplesheet <- samplesheet[samplesheet$cell_type %in% contrast,]
samplesheet <- samplesheet[samplesheet$source != "lung2",]
counts <- counts[,rownames(samplesheet)]
out <- rowSums(counts)
out <- names(out[out==0])
counts <- counts[!rownames(counts) %in% out,]
samplesheet <- data.frame("IDs" = rownames(samplesheet),
                          "celltype" = samplesheet$source,
                          "status" = samplesheet$cell_type,
                          "patient" = samplesheet$patient)

cancerAtlas_Regulatory_vs_helper_tcells <- list("counts" = counts, "samplesheet" = samplesheet)
saveRDS(object = cancerAtlas_Regulatory_vs_helper_tcells, file = "./cancerAtlas_Regulatory_vs_helper_tcells.rds")

