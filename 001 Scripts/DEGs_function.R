#' DEGs
#' @param data [[data.frame]] RNAseq data count matrix. Columns are cells and rows are the transcripts
#' @param samplesheet [[data.frame]] samplesheet with 4 columns. IDs, celltype, status and patient
#' @param cells [[vector]] Character vector of celltype names that need to be tested
#' @param contrast [[vector]] Character vector of the  two conditions that will be used for the contrast. First entry is reference
#' @return [[list]] 
DEGs <- function(data, samplesheet, cells, contrast, test="wilcox"){
  require(Seurat)  
  message("Preparing data...")
  S.obj <- CreateSeuratObject(counts = data)
  S.obj <- NormalizeData(S.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  S.obj@meta.data$celltype <- samplesheet$celltype
  S.obj@meta.data$status  <- samplesheet$status
  message("Starting tests...")
  res <- lapply(cells,function(cluster){
    message(cluster)
    tmp_subset <- subset(S.obj, subset = celltype  == cluster)    
    DE <- FindMarkers(tmp_subset, ident.1 = contrast[2], ident.2 = contrast[1], group.by = "status",logfc.threshold = 0,test.use = test)
    DE$fdr <- p.adjust(DE$p_val,method = "fdr")
    return(DE)
  })
  names(res) <- cells
  rm(data)
  rm(samplesheet)
  gc()
  return(res)
}