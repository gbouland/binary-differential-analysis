library(Seurat)
library(magrittr)
source("./DDs_function_parallel.R")
source("./DEGs_function_parallel.R")

## Example execution on complete AD dataset ##

set <- readRDS("./003 Data/AD.rds")#Load processed AD dataset
samplesheet <- set$samplesheet
counts <- set$counts
samplesheet$celltype <- "all"

## Execute in parallel
DDs <- DDsPar(data = counts,
              samplesheet = samplesheet,
              cells = "all",
              contrast = c("ct","AD"),
              cores = 12,
              chunks = 48,
              out = "./results/AD_res_DDs.rds")

DEGs <- DEGsPar(data = counts,
                samplesheet = samplesheet,
                cells = "all",
                contrast = c("ct","AD"),
                test="wilcox",
                cores=12,
                out = "./results/AD_res_wilcox.rds",
                norm=TRUE)

## Example execution on individual cell types AD dataset ##

set <- readRDS("./data/AD.rds") #Load processed AD dataset
samplesheet <- set$samplesheet
counts <- set$counts
cells <- c("oligo","astro","OPC","neuron","endo","mg")#specify cells types

## Execute in parallel
DDs <- DDsPar(data = counts,
              samplesheet = samplesheet,
              cells = cells,
              contrast = c("ct","AD"),
              cores = 12,
              chunks = 48,
              out = "./results/AD_res_DDs.rds")

DEGs <- DEGsPar(data = counts,
                samplesheet = samplesheet,
                cells = cells,
                contrast = c("ct","AD"),
                test="wilcox",
                cores=12,
                out = "./results/AD_res_wilcox.rds",
                norm=TRUE)