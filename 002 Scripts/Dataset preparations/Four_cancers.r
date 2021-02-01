library(Seurat)
library(magrittr)
library(Matrix)

meta <- rio::import("../GSE139555_all_metadata.txt")
tail(meta)

test <- readMM("./GSM4143680_SAM24360639-rb1.matrix.mtx.gz")
files <- list.files()
fileIDs <- unique(sapply(files,function(file){unlist(strsplit(x = file, split = "[.]"))[1]}) %>% unname())

allMTX <- lapply(fileIDs,function(fileID){
    mtx <- readMM(sprintf("./%s.matrix.mtx.gz",fileID))
    genes <- rio::import(sprintf("./%s.genes.tsv.gz",fileID),head=FALSE)
    cols <- rio::import(sprintf("./%s.barcodes.tsv.gz",fileID),head=FALSE)
    cols$V2 <- paste0(toupper(strsplit(x = fileID,split = "-")[[1]][2]),"_",cols$V1)
    colnames(mtx) <- cols$V2
    rownames(mtx) <- genes$V2
    return(mtx)
})

## Faster than do.call ##
rawCounts <- cbind(allMTX[[1]],allMTX[[2]],allMTX[[3]],allMTX[[4]],allMTX[[5]],allMTX[[6]],
              allMTX[[7]],allMTX[[8]],allMTX[[9]],allMTX[[10]],allMTX[[11]],allMTX[[12]],
              allMTX[[13]],allMTX[[14]],allMTX[[15]],allMTX[[16]],allMTX[[17]],allMTX[[18]],
              allMTX[[19]],allMTX[[20]],allMTX[[21]],allMTX[[22]],allMTX[[23]],allMTX[[24]],
              allMTX[[25]],allMTX[[26]],allMTX[[27]],allMTX[[28]],allMTX[[29]],allMTX[[30]],
              allMTX[[31]],allMTX[[32]])

head(meta)

## T-cells cluster identities## Source: https://www.nature.com/articles/s41586-020-2056-8.: Supplementary Table 2
tcell_clusters <- c(0,1,2,3,5,6,7,8,10,11,12,16,18,19,21,23,24,27)
meta <- meta[meta$ident %in% tcell_clusters,]
meta <- meta[meta$source %in% c("NAT","Tumor"),]


meta$celltype <- gsub('[[:digit:]]+', '', meta$patient)
samplesheet <- meta[,c("V1","celltype","source","patient")]
colnames(samplesheet) <- c("IDs","celltype","status","patient")

rawCounts <- rawCounts[,samplesheet$IDs]
out <- rowSums(rawCounts)
out <- names(out[out==0])
rawCounts <- rawCounts[!rownames(rawCounts) %in% out,]

fourCancers <- list("counts" = rawCounts, "samplesheet" = samplesheet)
saveRDS(object = fourCancers, file = "./fourCancers.rds")
