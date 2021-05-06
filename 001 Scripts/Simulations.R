library(muscat)
library(Seurat)
library(DDA)
library(SingleCellExperiment)
library(magrittr)

### wrapper function to simulate data and output as seurat object ###
generateData <- function(ref, n_similations, n_cells, n_genes, pct.de){
  known_simulations <- vector("list", n_similations)
  for(i in 1: n_similations){
    sim <- simData(ref, p_dd = c(1-pct.de,0,pct.de,0,0,0),
                   nc = n_cells, ng = n_genes, force = TRUE,nk = 1)
    counts <- sim@assays@data$counts
    gi <- metadata(sim)$gene_info
    samplesheet <- data.frame("colnames" = colnames(counts),
                              "group" = sim@colData@listData$group_id)
    rownames(counts) <- paste0(rownames(counts),"-",gi$category)
    obj <- CreateSeuratObject(counts = counts, project = "sim")
    obj[['contrast']] <- as.character(samplesheet$group)
    known_simulations[[i]] <- obj
  }
  return(known_simulations)
}

runTests <- function(seuratObject){
  ##DEA##
  seuratObject <- NormalizeData(seuratObject,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)
  DE_res <- FindMarkers(seuratObject,
                    ident.1 = "A",
                    ident.2 = "B",
                    group.by = "contrast",
                    logfc.threshold = 0,
                    test.use = "wilcox")
  DE_res$padjust <- p.adjust(DE_res$p_val, method = "fdr")
  ##DDA##
  
  seuratObject_DDA <- CreateSeuratObject(counts = seuratObject@assays$RNA@counts[rownames(DE_res),])
  seuratObject_DDA[['contrast']] <- seuratObject[['contrast']]
  DDA_res <- DDA(object = seuratObject_DDA,
             contrast.by = "contrast",
             ident.1 = "A",
             ident.2 = "B",
             p.adjust.method = "fdr")
  return(list("DEA" = DE_res,"DDA" = DDA_res$A_vs_B))
}

getPeformanceMetrics <- function(res){
  res$DE <- ifelse(grepl("-de", rownames(res)), 1, 0)
  res$sign <- ifelse(res$padjust <= 0.05,  1, 0)
  res$DE2 <- ifelse(grepl("-de", rownames(res)), "DE", "noDE")
  res$sign2 <- ifelse(res$padjust <= 0.05,  "Sign", "noSign")
  tableRes <- table(res$DE2, res$sign2)
  tp <- tableRes["DE", "Sign"]
  p <- sum(tableRes["DE", ])
  sens <- tp / p
  F1 <- MLmetrics::F1_Score(res$DE, res$sign, positive = "1")
  output <- data.frame(F1 = F1,
                       sens = sens,
                       "tp" = tableRes["DE", "Sign"],
                       "fp" = tableRes["noDE", "Sign"],
                       "tn" = tableRes["noDE", "noSign"],
                       "fn" = tableRes["DE", "noSign"])
  return(output)
}

### Provided reference ###
data(sce)
ref <- prepSim(sce, verbose = FALSE)

## Run simulations ##
sim_peform <- lapply(seq(from = 500, to = 2000, by = 500), function(n_cells){
  simulated <- generateData(ref = ref,
                            n_similations = 25,
                            n_cells = n_cells,
                            n_genes = 1000,
                            pct.de = 0.25)
  results <- lapply(simulated,runTests)
  metrics <- lapply(results, function(res){
    DDA_metrics <- getPeformanceMetrics(res$DDA)
    DEA_metrics <- getPeformanceMetrics(res$DEA)
    colnames(DDA_metrics) <- paste0("DDA_",colnames(DDA_metrics))
    colnames(DEA_metrics) <- paste0("DEA_",colnames(DEA_metrics))
    out <- cbind(DDA_metrics,DEA_metrics)
    return(out)
  }) %>% do.call(what = rbind)
  metrics$n_cells <- n_cells
  return(metrics)
}) %>% do.call(what = rbind)


## Load simulations results ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sim_peform <- readRDS("../003 Results/simulations_results.rds")

library(ggplot2)
##Accuracy
sel <- sim_peform[,c("DDA_F1","DEA_F1","n_cells")]
plotData <- lapply(unique(sel$n_cells),function(n){
  sumDDA <- summary(sel[sel$n_cells == n, "DDA_F1"])
  sumDEA <- summary(sel[sel$n_cells == n, "DEA_F1"])
  out <- list("DDA" = c("Test" = "DDA","median" = sumDDA[3],"Q1" = sumDDA[2],"Q3" = sumDDA[5],"n_cells" = n),
              "DEA" = c("Test" = "DEA","median" = sumDEA[3],"Q1" = sumDEA[2],"Q3" = sumDEA[5],"n_cells" = n))
  out <- do.call(what = rbind,out)
  return(out)
}) %>% do.call(what = rbind) %>% as.data.frame()
colnames(plotData) <- c("Test","median","Q1","Q3","n_cells")
plotData$Test <- as.factor(plotData$Test)
plotData$median <- as.numeric(plotData$median)
plotData$Q1<- as.numeric(plotData$Q1)
plotData$Q3<- as.numeric(plotData$Q3)
plotData$n_cells<- factor(plotData$n_cells, levels = c("500","1000","1500","2000"))
acc <- ggplot(plotData, aes(Test,median,fill = Test)) + geom_bar(stat = "identity") + geom_errorbar( aes(x=Test, ymin=Q1, ymax=Q3), width=0.4, alpha=0.9, size=1)+
  facet_grid(~n_cells) + ylim(0,1) + scale_fill_manual(values=c("#56B4E9", "#009E73")) + theme_minimal() + ylab("F1 score") + theme(legend.position = "none")


#Sensitivity
sel <- sim_peform[,c("DDA_sens","DEA_sens","n_cells")]
colnames(sel) <- c("DDA","DEA","n_cells")
plotData <- lapply(unique(sel$n_cells),function(n){
  sumDDA <- summary(sel[sel$n_cells == n, "DDA"])
  sumDEA <- summary(sel[sel$n_cells == n, "DEA"])
  out <- list("DDA" = c("Test" = "DDA","median" = sumDDA[3],"Q1" = sumDDA[2],"Q3" = sumDDA[5],"n_cells" = n),
              "DEA" = c("Test" = "DEA","median" = sumDEA[3],"Q1" = sumDEA[2],"Q3" = sumDEA[5],"n_cells" = n))
  out <- do.call(what = rbind,out)
  return(out)
}) %>% do.call(what = rbind) %>% as.data.frame()
colnames(plotData) <- c("Test","median","Q1","Q3","n_cells")
plotData$Test <- as.factor(plotData$Test)
plotData$median <- as.numeric(plotData$median)
plotData$Q1<- as.numeric(plotData$Q1)
plotData$Q3<- as.numeric(plotData$Q3)
plotData$n_cells<- factor(plotData$n_cells, levels = c("500","1000","1500","2000"))
sens <- ggplot(plotData, aes(Test,median,fill = Test)) + geom_bar(stat = "identity") + geom_errorbar( aes(x=Test, ymin=Q1, ymax=Q3), width=0.4, alpha=0.9, size=1)+
  facet_grid(~n_cells) + ylim(0,1) + scale_fill_manual(values=c("#56B4E9", "#009E73")) + theme_minimal() + ylab("TPR") + theme(legend.position = "none")

##FPR
sel <- apply(sim_peform,1,function(row){
  DDA <- row["DDA_fp"] / (row["DDA_fp"] + row["DDA_tn"])
  DEA <- row["DEA_fp"] / (row["DEA_fp"] + row["DEA_tn"])
  return(c(DDA,DEA,row["n_cells"]))
}) %>% t() %>% as.data.frame()
colnames(sel) <-c("DDA","DEA","n_cells")
plotData <- lapply(unique(sel$n_cells),function(n){
  sumDDA <- summary(sel[sel$n_cells == n, "DDA"])
  sumDEA <- summary(sel[sel$n_cells == n, "DEA"])
  out <- list("DDA" = c("Test" = "DDA","median" = sumDDA[3],"Q1" = sumDDA[2],"Q3" = sumDDA[5],"n_cells" = n),
              "DEA" = c("Test" = "DEA","median" = sumDEA[3],"Q1" = sumDEA[2],"Q3" = sumDEA[5],"n_cells" = n))
  out <- do.call(what = rbind,out)
  return(out)
}) %>% do.call(what = rbind) %>% as.data.frame()
colnames(plotData) <- c("Test","median","Q1","Q3","n_cells")
plotData$Test <- as.factor(plotData$Test)
plotData$median <- as.numeric(plotData$median)
plotData$Q1<- as.numeric(plotData$Q1)
plotData$Q3<- as.numeric(plotData$Q3)
plotData$n_cells<- factor(plotData$n_cells, levels = c("500","1000","1500","2000"))
fpr <- ggplot(plotData, aes(Test,median,fill = Test)) + geom_bar(stat = "identity") + geom_errorbar( aes(x=Test, ymin=Q1, ymax=Q3), width=0.4, alpha=0.9, size=1)+
  facet_grid(~n_cells)  + scale_fill_manual(values=c("#56B4E9", "#009E73")) + theme_minimal() + ylab("FPR") + theme(legend.position = "none") #+ ylim(0,0.5)


##PPV
sel <- apply(sim_peform,1,function(row){
  DDA <- row["DDA_tp"] / (row["DDA_tp"] + row["DDA_fp"])
  DEA <- row["DEA_tp"] / (row["DEA_tp"] + row["DEA_fp"])
  return(c(DDA,DEA,row["n_cells"]))
}) %>% t() %>% as.data.frame()
colnames(sel) <-c("DDA","DEA","n_cells")
plotData <- lapply(unique(sel$n_cells),function(n){
  sumDDA <- summary(sel[sel$n_cells == n, "DDA"])
  sumDEA <- summary(sel[sel$n_cells == n, "DEA"])
  out <- list("DDA" = c("Test" = "DDA","median" = sumDDA[3],"Q1" = sumDDA[2],"Q3" = sumDDA[5],"n_cells" = n),
              "DEA" = c("Test" = "DEA","median" = sumDEA[3],"Q1" = sumDEA[2],"Q3" = sumDEA[5],"n_cells" = n))
  out <- do.call(what = rbind,out)
  return(out)
}) %>% do.call(what = rbind) %>% as.data.frame()
colnames(plotData) <- c("Test","median","Q1","Q3","n_cells")
plotData$Test <- as.factor(plotData$Test)
plotData$median <- as.numeric(plotData$median)
plotData$Q1<- as.numeric(plotData$Q1)
plotData$Q3<- as.numeric(plotData$Q3)
plotData$n_cells<- factor(plotData$n_cells, levels = c("500","1000","1500","2000"))
ppv <- ggplot(plotData, aes(Test,median,fill = Test)) + geom_bar(stat = "identity") + geom_errorbar( aes(x=Test, ymin=Q1, ymax=Q3), width=0.4, alpha=0.9, size=1)+
  facet_grid(~n_cells)  + scale_fill_manual(values=c("#56B4E9", "#009E73")) + theme_minimal() + ylab("PPV") + ylim(0,1) + theme(legend.position = "none")






