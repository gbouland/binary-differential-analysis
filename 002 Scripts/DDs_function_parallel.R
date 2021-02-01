#' DDsPar
#' @param data [[data.frame]] RNAseq data count matrix. Columns are cells and rows are the transcripts
#' @param samplesheet [[data.frame]] samplesheet with 4 columns. IDs, celltype, status and patient
#' @param cells [[vector]] Character vector of celltype names that need to be tested
#' @param contrast [[vector]] Character vector of the  two conditions that will be used for the contrast. First entry is reference
#' @return [[list]] 
DDsPar <- function(data, samplesheet, cells, contrast, cores, chunks, out){    
  require(doParallel)    
  geneList <- split(sample(rownames(data)), sort(1:length(rownames(data))%%chunks))    
  subCounts <- vector("list", length = chunks)
  start_time <- Sys.time()
  for(i in 1:chunks){
    message(i)
    subCounts[[i]] <- as.matrix(data[geneList[[i]],])
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  rm(data)
  rm(geneList)
  gc()
  cl <- makeCluster(cores, type="PSOCK",outfile="./log.txt")
  registerDoParallel(cl)
  res <- foreach(sub=subCounts) %dopar% {
    source("./DDs_function.R")
    DDs(data = sub, samplesheet = samplesheet, cells = cells, contrast = contrast)
  }
  stopCluster(cl)
  res <- do.call(what = 'c',res) 
  backTogether <- lapply(unique(names(res)),function(name){
    tmp.togh <- do.call(what = rbind, res[names(res) == name])
    rownames(tmp.togh) <- make.unique(gsub(paste0(name,"."),"",rownames(tmp.togh)))
    tmp.togh$fdr <- p.adjust(tmp.togh$P,method = "fdr")
    tmp.togh <- tmp.togh[order(tmp.togh$P),]
    return(tmp.togh)
  })
  names(backTogether) <- unique(names(res))
  #backTogether
  saveRDS(backTogether, file = out)
}