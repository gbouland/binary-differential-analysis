generateDataIndepth <- function(data,summ,set){
  all_DD <- lapply(data,function(x){x$DD}) %>% do.call(what= "rbind")
  all_will <- lapply(data,function(x){x$wilcox}) %>% do.call(what= "rbind")
  common <- intersect(rownames(all_DD),rownames(all_will))
  all_DD <- all_DD[common,]
  all_will <- all_will[common,]
  x <- data.frame("gene" = common,
                  "P_DD" = all_DD$P,
                  "P_wil" = all_will$p_val,
                  "P_delta" = abs(-log10(all_DD$P) - -log10(all_will$p_val)),
                  "DD_sign" = ifelse(all_DD$fdr<=0.05,1,0),
                  "wil_sign" = ifelse(all_will$fdr<=0.05,1,0))
  x$dataset <- strsplit(x$gene,split = "[.]") %>% sapply(FUN=function(x)x[1]) %>% unname()
  x$gene <- strsplit(x$gene,split = "[.]") %>% sapply(FUN=function(x)x[2]) %>% unname()
  x$P_dif <- ifelse(x$P_delta >= 5, "Big difference", "small difference")
  x$set <- ifelse(x$DD_sign ==1 & x$wil_sign ==1,"common","not sign")
  x$set <- ifelse(x$DD_sign ==1 & x$wil_sign ==0,"DD only",x$set)
  x$set <- ifelse(x$DD_sign ==0 & x$wil_sign ==1,"Wil only",x$set)
  x$P_dif <- ifelse(x$set == "common", NA, x$P_dif)
  toExplore <- x[x$dataset == set,]
  allStats <- summ[[set]]
  allStats$DD_sign <- toExplore[match(rownames(allStats),toExplore$gene),"DD_sign"]
  allStats$wil_sign <- toExplore[match(rownames(allStats),toExplore$gene),"wil_sign"]
  allStats$P_DD <- toExplore[match(rownames(allStats),toExplore$gene),"P_DD"]
  allStats$P_wil <- toExplore[match(rownames(allStats),toExplore$gene),"P_wil"]
  allStats$P_delta <- toExplore[match(rownames(allStats),toExplore$gene),"P_delta"]
  allStats <- allStats[!is.na(allStats$DD_sign) | !is.na(allStats$wil_sign), ]
  allStats$set <- ifelse(allStats$DD_sign ==1 & allStats$wil_sign ==1,"common","not sign")
  allStats$set <- ifelse(allStats$DD_sign ==1 & allStats$wil_sign ==0,"DD only",allStats$set)
  allStats$set <- ifelse(allStats$DD_sign ==0 & allStats$wil_sign ==1,"Wil only",allStats$set)
  allStats$P_dif <- ifelse(allStats$P_delta >= 5, "Big difference", "small difference")
  allStats$gene <- rownames(allStats)
  return(allStats)
}