Pancreas <- rio::import("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/Pancreas/pancreas_refseq_rpkms_counts_3514sc.txt")
cols <- colnames(Pancreas)[1:3516]
Pancreas <- Pancreas[,c(1:2,3517:7030)]
colnames(Pancreas) <- cols
Pancreas$Accession <- NULL
geneNames <- Pancreas$Gene
Pancreas$Gene <- NULL
samp <- rio::import("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/Pancreas/E-MTAB-5061.sdrf.txt")
Pancreas <- t(Pancreas)
colnames(Pancreas) <- geneNames
samp <- samp[samp$`Characteristics[single cell well quality]` == "OK",]
Pancreas <- Pancreas[samp$`Source Name`,]
Pancreas <- t(Pancreas)
samp <- samp[,c("Source Name","Characteristics[cell type]","Characteristics[disease]","Characteristics[individual]")]
colnames(samp) <- c("IDs","celltype","status","patient")
T2Dset <- list("counts" = Pancreas, "samplesheet" = samp)
saveRDS(object = T2Dset, file = "C:/Users/gabouland/Documents/004 PhD/008scBinaryAnalysis/000 Preprocess/T2Dset.rds")
