library(Matrix)
MDD <- readMM("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/MDD/GSE144136_GeneBarcodeMatrix_Annotated.mtx")
geneNames <- rio::import("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/MDD/GSE144136_GeneNames.csv")
geneNames<- geneNames[2:nrow(geneNames),]
geneNames$V1 <- NULL
rownames(MDD) <- geneNames$V2
load("C:/Users/gabouland/Documents/004 PhD/000Data/SingleCell/MDD/samplesheet.RData")
colnames(MDD) <- samplesheet$ColN
samplesheet <- samplesheet[,c("ColN","cellType","Group","SubjectID")]
colnames(samplesheet) <- c("IDs","celltype","status","patient")
samplesheet$celltype <- ifelse(grepl("Oligos",samplesheet$celltype),"Oligos",samplesheet$celltype)
samplesheet$celltype <- ifelse(grepl("Astros",samplesheet$celltype),"Astros",samplesheet$celltype)
samplesheet$celltype <- ifelse(grepl("Ex_",samplesheet$celltype),"Ex",samplesheet$celltype)
samplesheet$celltype <- ifelse(grepl("Inhib_",samplesheet$celltype),"Inhib",samplesheet$celltype)
samplesheet$celltype <- ifelse(grepl("OPCs_",samplesheet$celltype),"OPCs",samplesheet$celltype)
cells.MDD <- c("Astros","Endo","Micro/Macro","Oligos","OPCs")
samplesheet <- samplesheet[samplesheet$celltype %in% cells.MDD,]
MDD <- MDD[,samplesheet$IDs]
MDD_set <- list("counts" = MDD, "samplesheet" = samplesheet)
saveRDS(object = MDD_set, file = "C:/Users/gabouland/Documents/004 PhD/008scBinaryAnalysis/000 Preprocess/MDD_set.rds")
