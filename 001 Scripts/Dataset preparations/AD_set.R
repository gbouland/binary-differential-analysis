rawCounts <- rio::import("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/AD_set/scRNA_rawCounts.tsv")
samplesheet <- rio::import("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/AD_set/scRNA_metadata.tsv")
rownames(rawCounts) <- rawCounts$geneName
rawCounts$geneName <- NULL
samplesheet <- samplesheet[,c("sampleID","cellType","batchCond","patient")]
colnames(samplesheet) <- c("IDs","celltype","status","patient")
rawCounts <- rawCounts[,samplesheet$IDs]


ADset <- list("counts" = rawCounts, "samplesheet" = samplesheet)

saveRDS(object = ADset, file = "C:/Users/gabouland/Documents/004 PhD/008scBinaryAnalysis/000 Preprocess/ADset.rds")



head(colnames(rawCounts))



samplesheet
