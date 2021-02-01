#' DEGsPar
#' @param data [[data.frame]] RNAseq data count matrix. Columns are cells and rows are the transcripts
#' @param samplesheet [[data.frame]] samplesheet with 4 columns. IDs, celltype, status and patient
#' @param cells [[vector]] Character vector of celltype names that need to be tested
#' @param contrast [[vector]] Character vector of the  two conditions that will be used for the contrast. First entry is reference
#' @return [[list]] 
DEGsPar <- function(data, samplesheet, cells, contrast, test="wilcox", cores, out, norm){    
  require(doParallel) 
  require(Seurat)
  require(magrittr)
  S.obj <- CreateSeuratObject(counts = data)
  if(norm){
    S.obj <- NormalizeData(S.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  }    
  geneList <- split(sample(rownames(S.obj)), sort(1:length(rownames(S.obj))%%cores))
  start_time <- Sys.time()
  sub.Objs <- vector("list", length = cores)
  for(i in 1:cores){
    message(i)
    sub.Objs[[i]] <- CreateSeuratObject(counts = S.obj@assays$RNA@data[geneList[[i]],])
  }
  end_time <- Sys.time()
  print(end_time - start_time)  
  rm(S.obj)
  rm(data)
  rm(geneList)
  gc()    
  cl <- makeCluster(cores, type="PSOCK",outfile="/home/nfs/gabouland/log.txt")
  registerDoParallel(cl)
  message("start parallel...")
  start_time <- Sys.time()
  res <- foreach(tmp.obj=sub.Objs) %dopar% {     
    require(Seurat)
    require(magrittr)
    message("In parallel...")
    message(dim(tmp.obj@assays$RNA@counts))
    tmp.obj@meta.data$celltype <- samplesheet$celltype
    tmp.obj@meta.data$status  <- samplesheet$status
    cell_res <- lapply(cells,function(cluster){
      message(cluster)
      tmp_subset <- subset(tmp.obj, subset = celltype  == cluster)    
      DE <- FindMarkers(tmp_subset, ident.1 = contrast[2], ident.2 = contrast[1], group.by = "status",logfc.threshold = 0, test.use = test)
      #DE$fdr <- p.adjust(DE$p_val,method = "fdr")
      return(DE)
    })
    names(cell_res) <- cells
    message("finished..")
    cell_res
  } %>% do.call(what = 'c')    
  stopCluster(cl)
  backTogether <- lapply(unique(names(res)),function(name){
    tmp.togh <- do.call(what = rbind, res[names(res) == name])
    rownames(tmp.togh) <- make.unique(gsub(paste0(name,"."),"",rownames(tmp.togh)))
    tmp.togh$fdr <- p.adjust(tmp.togh$p_val,method = "fdr")
    tmp.togh <- tmp.togh[order(tmp.togh$p_val),]
    return(tmp.togh)
  })
  names(backTogether) <- unique(names(res))
  end_time <- Sys.time()
  print(end_time - start_time)
  saveRDS(backTogether, file = out)
  #return(backTogether)
}