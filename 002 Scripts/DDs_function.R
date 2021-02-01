#' DDs
#' @param data [[data.frame]] RNAseq data count matrix. Columns are cells and rows are the transcripts
#' @param samplesheet [[data.frame]] samplesheet with 4 columns. IDs, celltype, status and patient
#' @param cells [[vector]] Character vector of celltype names that need to be tested
#' @param contrast [[vector]] Character vector of the  two conditions that will be used for the contrast. First entry is reference
#' @return [[list]] 
DDs <- function(data, samplesheet, cells, contrast){
  DDres <- lapply(cells,function(cell){
    message(cell)    
    cell.samplesheet <- samplesheet[samplesheet$celltype == cell,]
    message("binarizing...")
    cell.Counts <- as.matrix(data[,cell.samplesheet$IDs])
    cell.Counts <- cell.Counts >= 1    
    results <- matrix(nrow = nrow(cell.Counts),ncol = 8)
    total <- nrow(cell.Counts)
    message("Starting tests...")
    for(i in 1:total){      
      if(i %% 1000 == 0){message(sprintf("%s of %s",i,total))}      
      tmp <- data.frame("Gene" = cell.Counts[i,], "status" = cell.samplesheet$status)
      tmp$status <- factor(tmp$status,levels = contrast)
      tab <- table(tmp[,1],tmp[,2])
      N1 <- sum(tab[,1])
      N2 <- sum(tab[,2])
      pct1 <- tab[1,1] / N1 
      pct2 <- tab[1,2] / N2
      if(pct1 < 0.999 | pct2 < 0.999){ 
        res <- glm("status ~ Gene",family = "binomial",data = tmp)
        sum <- summary(res)
        newline <- sum$coefficients[2,]
        newline <- c(newline,c(round(pct1,3),round(pct2,2),N1,N2))               
      }
      else{     
        newline <- c(NA,NA,NA,NA,NA,NA,NA,NA)        
      }      
      results[i,] <- newline      
    }
    results <- data.frame(results)
    rownames(results) <- make.unique(rownames(cell.Counts))
    colnames(results) <- c("Estimate","Std.Error","Z","P","pct1","pct2","N1","N2")
    results <- na.omit(results)
    results$fdr <- p.adjust(results$P,method="fdr")
    results <- results[order(results$fdr),]
    return(results)
  })
  names(DDres) <- cells
  return(DDres)
}