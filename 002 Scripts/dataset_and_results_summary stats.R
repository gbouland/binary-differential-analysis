getDataStats <- function(counts,samplesheet, contrast){
  control_ids <- as.character(samplesheet[samplesheet$status == contrast[1],"IDs"])
  case_ids <- as.character(samplesheet[samplesheet$status == contrast[2],"IDs"])
  zeroMeasurements_controls <- (as.double(sum(counts[,control_ids] == 0))) / (as.double(nrow(counts[,control_ids])) * as.double(ncol(counts[,control_ids])))
  zeroMeasurements_cases <- (as.double(sum(counts[,case_ids] == 0))) / (as.double(nrow(counts[,case_ids])) * as.double(ncol(counts[,case_ids])))
  
  NInd_controls <- length(unique(samplesheet[samplesheet$IDs %in% control_ids,"patient"]))
  NInd_cases <- length(unique(samplesheet[samplesheet$IDs %in% case_ids,"patient"]))
  
  Nsubpopu <- length(unique(samplesheet$celltype))
  Ngenes <- nrow(counts)
  Ncells <- ncol(counts)
  
  data.frame("zeroPct_ct" = zeroMeasurements_controls,
             "zeroPct_case" = zeroMeasurements_cases,
             "N_individuals_ct" = NInd_controls,
             "N_individuals_case" = NInd_cases,
             "No.cells_ct" = length(control_ids),
             "No.cells_case" = length(case_ids),
             "N_subpopulations" = Nsubpopu,
             "N_genes" = Ngenes,
             "N_cells" = Ncells)
}

getResultsStats <- function(res){
  DDres <- res$DD
  wilcoxres <- res$wilcox
  DDres_sign <- DDres[DDres$fdr <= 0.05,]
  wilcoxres_sign <- wilcoxres[wilcoxres$fdr <= 0.05,]
  
  DDres_genes <- rownames(DDres)
  DDres_sign_genes <- rownames(DDres_sign)
  
  wilcoxres_genes <- rownames(wilcoxres)
  wilcoxres_sign_genes <- rownames(wilcoxres_sign)
  
  N_total_tested <- length(unique(c(DDres_genes,wilcoxres_genes)))    
  
  N_DD_unique_tested <- length(DDres_genes[!DDres_genes %in% wilcoxres_genes])
  N_wilcox_unique_tested <- length(wilcoxres_genes[!wilcoxres_genes %in% DDres_genes])
  
  N_total_sign <- length(unique(c(DDres_sign_genes,wilcoxres_sign_genes)))
  N_DD_sign <- length(DDres_sign_genes)
  N_wilcox_sign <- length(wilcoxres_sign_genes)
  N_common_sign <- length(intersect(DDres_sign_genes,wilcoxres_sign_genes))
  
  N_DD_unique_sign <- length(DDres_sign_genes[!DDres_sign_genes %in% wilcoxres_sign_genes])
  N_wilcox_unique_sign <- length(wilcoxres_sign_genes[!wilcoxres_sign_genes %in% DDres_sign_genes])
  
  pctAgree <- N_common_sign / N_total_sign    
  
  cor <-cor.test(DDres$Estimate, wilcoxres[rownames(DDres),"avg_logFC"])$estimate 
  
  data.frame("N_total_tested" = N_total_tested,
             "N_total_sign" = N_total_sign,
             "N_common_sign" = N_common_sign,
             "N_DD_sign" = N_DD_sign,
             "N_wilcox_sign" = N_wilcox_sign,
             "N_DD_unique_sign" = N_DD_unique_sign,
             "N_wilcox_unique_sign" = N_wilcox_unique_sign,
             "pct agree" = pctAgree,
             "cor" = cor)
}

## Code to get the dataset statistics and results statitics / comparison ##
##Example based on AD dataset##

library(magrittr)
nothreshDat <- readRDS("./results/DDA_DEA_results_all_datasets.rds")
set <- readRDS("./002 Preprocessed/AD.rds")
samplesheet <- set$samplesheet
counts <- set$counts
counts <- counts[,samplesheet$IDs]
AD <- cbind(getDataStats(counts, samplesheet, c("ct","AD")),getResultsStats(nothreshDat$AD))
