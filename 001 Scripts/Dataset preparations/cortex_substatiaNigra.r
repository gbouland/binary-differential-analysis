library(Seurat)
library(magrittr)
samplesheet <- rio::import("/tudelft.net/staff-umbrella/pQTL/000scBinary/000 Data/011cortex_substatiaNigra/samplesheet.csv")
setwd("/tudelft.net/staff-umbrella/pQTL/000scBinary/000 Data/011cortex_substatiaNigra/")

dirs <- list.dirs()
dirs <- dirs[grep("GSM",dirs)]
dirs <- gsub("./","",dirs)


all_mtx <- lapply(dirs,function(dir){
    print(dir)
    batch <- unlist(strsplit(split = "_", dir))[4]
    expression_matrix <- Read10X(data.dir = sprintf("/tudelft.net/staff-umbrella/pQTL/000scBinary/000 Data/011cortex_substatiaNigra/%s",dir))
    cols <- colnames(expression_matrix)
    cols <- paste0(batch,"_",cols)
    cols <- gsub("-1","",cols)    
    colnames(expression_matrix) <- cols
    return(expression_matrix)
})
rawCounts <- cbind(all_mtx[[1]],all_mtx[[2]],all_mtx[[3]],all_mtx[[4]],all_mtx[[5]],all_mtx[[6]],
              all_mtx[[7]],all_mtx[[8]],all_mtx[[9]],all_mtx[[10]],all_mtx[[11]],all_mtx[[12]])

rawCounts <- rawCounts[,samplesheet$IDs]
sums <- as.matrix(rawCounts)
sums <- rowSums(sums)
exclude <- names(sums[sums==0])
rawCounts <- rawCounts[!rownames(rawCounts) %in% exclude,]


table(samplesheet$Brain_Region,samplesheet$Level_1_cell_type)


### CORTEX ###
samplesheet_cortex <- samplesheet[samplesheet$Brain_Region == "Cortex",]
samplesheet_cortex <- samplesheet_cortex[samplesheet_cortex$Level_1_cell_type %in% c("Excitatory neurons","Inhibitory neurons","ODC","OPC"),]
samplesheet_cortex$celltype <- ifelse(grepl("neurons",samplesheet_cortex$Level_1_cell_type),"neurons",NA)
samplesheet_cortex$celltype <- ifelse(samplesheet_cortex$Level_1_cell_type %in% c("OPC","ODC"),"oligo",samplesheet_cortex$celltype)
samplesheet_cortex <- samplesheet_cortex[,c("IDs","celltype","Level_1_cell_type","Library")]
colnames(samplesheet_cortex) <- c("IDs","celltype","status","patient")
counts_cortex <- rawCounts[,samplesheet_cortex$IDs]
cortex_set <- list("counts" = counts_cortex, "samplesheet" = samplesheet_cortex)
saveRDS(object = cortex_set, file = "/tudelft.net/staff-umbrella/pQTL/000scBinary/002 Preprocessed/cortex_set.rds")


### substantia nigra ###
samplesheet_subs <- samplesheet[samplesheet$Brain_Region == "Substantia nigra",]
samplesheet_subs <- samplesheet_subs[samplesheet_subs$Level_1_cell_type %in% c("ODC","OPC"),]
samplesheet_subs$celltype <- "oligo"
samplesheet_subs <- samplesheet_subs[,c("IDs","celltype","Level_1_cell_type","Library")]
colnames(samplesheet_subs) <- c("IDs","celltype","status","patient")
counts_subs <- rawCounts[,samplesheet_subs$IDs]

substantiaNigra_set <- list("counts" = counts_subs, "samplesheet" = samplesheet_subs)
saveRDS(object = substantiaNigra_set, file = "/tudelft.net/staff-umbrella/pQTL/000scBinary/002 Preprocessed/substantiaNigra_set.rds")
