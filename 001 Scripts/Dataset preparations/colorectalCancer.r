library(magrittr)
counts <- rio::import("/tudelft.net/staff-umbrella/pQTL/000scBinary/000 Data/010colorectalCancer/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt")
samplesheet <- rio::import("/tudelft.net/staff-umbrella/pQTL/000scBinary/000 Data/010colorectalCancer/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz")

rownames(counts) <- counts$Index
counts$Index <- NULL

samplesheet <- samplesheet[,c("Index","Cell_type","Class","Patient")]
colnames(samplesheet) <- c("IDs","cellype","status","patient")
selection <- table(samplesheet$cellype,samplesheet$status) %>% as.matrix()
selection <- selection[selection[,1]>5  & selection[,2]>5,] %>% rownames()
samplesheet <- samplesheet[samplesheet$cellype %in% selection,]
counts <- counts[,samplesheet$IDs]
sums <- rowSums(counts)
exclude <- names(sums[sums == 0])
counts <- counts[!rownames(counts) %in% exclude,]


colorectalCancer_set <- list("counts" = counts, "samplesheet" = samplesheet)
saveRDS(object = colorectalCancer_set, file = "/tudelft.net/staff-umbrella/pQTL/000scBinary/002 Preprocessed/colorectalCancer.rds")


