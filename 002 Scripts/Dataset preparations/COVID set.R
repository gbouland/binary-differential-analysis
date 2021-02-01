covid <- readRDS("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/Covid19/blish_covid.seu.rds")
samplesheet <- covid@meta.data
samplesheet <- samplesheet[,c("orig.ident","nCount_RNA","nFeature_RNA","cell.type.fine","cell.type.coarse","cell.type","Status","Donor.full","Donor","Sex","Ventilated")]
rawCounts <- covid@assays$RNA@counts

samplesheet <- samplesheet[,c("cell.type.coarse","Status","Donor")]
samplesheet$IDs <- rownames(samplesheet)
colnames(samplesheet) <- c("celltype","status","patient","IDs")
COVID <- list("counts" = rawCounts, "samplesheet" = samplesheet)
saveRDS(object = COVID, file = "C:/Users/gabouland/Documents/004 PhD/008scBinaryAnalysis/000 Preprocess/COVID.rds")
  