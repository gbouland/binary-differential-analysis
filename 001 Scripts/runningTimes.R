library(muscat)
library(Seurat)
library(DDA)
library(SingleCellExperiment)
library(reshape2)


### Provided reference ###
data(sce)
ref <- prepSim(sce, verbose = FALSE)


n_genes <- 500
n_cells <- 10000

sim <- simData(ref,
               p_dd = c(0.75, 0, 0.25, 0, 0, 0),
               nc = n_cells,
               ng = n_genes,
               force = TRUE,
               nk = 1)

counts <- sim@assays@data$counts
gi <- metadata(sim)$gene_info
samplesheet <- data.frame("colnames" = colnames(counts),
                          "group" = sim@colData@listData$group_id)
rownames(counts) <- paste0(rownames(counts), "-", gi$category)


time_results <- matrix(nrow = 10, ncol = 10)
seurat_methods <- c("wilcox",
                    "bimod",
                    "roc",
                    "t",
                    "negbinom",
                    "poisson",
                    "LR",
                    "MAST",
                    "DESeq2")
colnames(time_results) <- c("DDA",seurat_methods)



for(h in 1:10){
  i <- seq(from = 1000, to = 10000, by = 1000)[h]
  sub <- counts[,sample(1:10000,size = i)]
  sub_samplesheet <- samplesheet[match(colnames(sub),samplesheet$colnames),]
  obj <- CreateSeuratObject(counts = sub, project = "sim")
  obj[['contrast']] <- as.character(sub_samplesheet$group)
  
  time_results[h,"DDA"] <- system.time({
    results <- DDA(object = obj,
                   contrast.by = "contrast",
                   ident.1 = "A",
                   ident.2 = "B",
                   p.adjust.method = "fdr",
                   min.pct = 0,
                   n.cores = 1
                   )
    })[3]
  
  obj <- NormalizeData(obj,
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)
  
  for(method in seurat_methods){
    message(method)
    time_results[h,method] <- system.time({
      DE <- FindMarkers(obj,
                        ident.1 = "A",
                        ident.2 = "B",
                        group.by = "contrast",
                        logfc.threshold = 0,
                        test.use = method,
                        min.pct = 0)
      })[3]
  }
  print(time_results)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
time_results <- readRDS("../003 Results/timing_res.rds")
rownames(time_results) <- paste0(seq(from = 1000, to = 10000, by = 1000))
time_results_df <- melt(time_results)
ggplot(time_results_df,aes(as.factor(Var1), value, group = Var2, col = Var2)) +
  geom_point() + geom_line(linetype = "dashed") + theme_minimal() + ylab("Run time (s)") +
  xlab("No. of cells") + labs(col = "Test")













