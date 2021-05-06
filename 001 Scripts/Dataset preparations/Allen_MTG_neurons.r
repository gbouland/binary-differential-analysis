
setwd("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/Allen_mtg")


exons <- rio::import("human_MTG_2018-06-14_exon-matrix.csv")
introns <- rio::import("human_MTG_2018-06-14_intron-matrix.csv")
both <- exons + introns
rm(introns);rm(exons)


genes <- rio::import("human_MTG_2018-06-14_genes-rows.csv")
rownames(both) <- genes$gene
both$V1 <- NULL
meta <- rio::import("human_MTG_2018-06-14_samples-columns.csv")
meta <- meta[,c("sample_name","brain_subregion","class","donor")]
colnames(meta) <- c("IDs","celltype","status","patient")
meta <- meta[meta$status %in% c("GABAergic","Glutamatergic"),]

both <- both[,meta$IDs]
out <- rowSums(both)
out <- names(out[out==0])
both <- both[!rownames(both) %in% out,]
MTG_neuron <- list("counts" = both, "samplesheet" = meta)
saveRDS(object = MTG_neuron, file = "C:/Users/gabouland/Documents/004 PhD/008scBinaryAnalysis/000 Preprocess/MTG_neuron.rds")


