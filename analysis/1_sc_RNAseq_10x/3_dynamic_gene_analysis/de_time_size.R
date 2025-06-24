### DE for aggressiveness in each cluster
rm(list = ls())
setwd("/singerlab/linglin/b16/")

library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(edgeR)
library(Matrix)

#### configuration ####
today <- "2020-12-30"
prep_data_date <- "2020-08-06"

dir_out <- paste0("2_pipeline/TSA/DE_time_size/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)


#### differential expression ####
cell_type_vec <- c("T", "nT")
for (cell_type in cell_type_vec) {
  ## load data
  so <- readRDS(paste0("2_pipeline/prep_data/so_", cell_type, "_", prep_data_date, ".rds"))
  
  ## compute pseudobulk count matrix
  so$sample_cluster <- paste(so$sample_num, so$louvain, sep = "_")
  sample_cluster_vec <- unique(so$sample_cluster) %>% sort()
  psb_mtx <- sapply(sample_cluster_vec, function(x) {
    idx <- which(so$sample_cluster == x)
    if (length(idx) == 1) return(so[["RNA"]]@counts[,idx]) ## only one cell is in the cluster
    so[["RNA"]]@counts[,idx] %>% Matrix::rowSums(.)
  })
  
  ## create pseudobulk meta data
  psb_meta <- strsplit(sample_cluster_vec, split = "_") %>% Reduce(rbind, .) %>% data.frame()
  rownames(psb_meta) <- sample_cluster_vec
  colnames(psb_meta) <- c("sample_num", "cluster")
  psb_meta$mouse_id <- so$mouse_id[match(psb_meta$sample_num, so$sample_num)]
  psb_meta$time <- so$time_point[match(psb_meta$sample_num, so$sample_num)] %>% as.numeric()
  psb_meta$size <- so$tumor_size[match(psb_meta$sample_num, so$sample_num)] %>% as.numeric()
  
  
  ## edgeR in each cluster
  cluster_vec <- unique(so$louvain) %>% sort()
  res_time <- list()
  res_size <- list()
  for (cluster in cluster_vec) {
    # if (cluster == 12 & cell_type == "non_T_cell") next() ## need to uncomment for factor time becasue not enough samples at each time point (9 11 12 13 16) = (1  0  1  1  3)
    cat("Running Pseudo-bulk EdgeR for", cell_type, "Cluster", cluster, "...\n")
    psb_meta_sub <- psb_meta[psb_meta$cluster == cluster,] %>% droplevels()
    
    if (nrow(psb_meta_sub) < 2) {
      res[[paste0("Cluster_", cluster)]] <- NULL
      next()
    }
    
    dge <- DGEList(counts = psb_mtx[,psb_meta$cluster == cluster])
    dge <- calcNormFactors(dge, method="TMM")
    
    ## time numeric
    design <- model.matrix(~ time, data = psb_meta_sub)
    dge <- estimateDisp(dge, design = design)
    fit <- glmFit(dge, design = design)
    lrt <- glmLRT(fit, coef = 2)
    res_time[[paste0("Cluster_", cluster)]] <- topTags(lrt, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)$table
    
    ## size numeric
    design <- model.matrix(~ size, data = psb_meta_sub)
    dge <- estimateDisp(dge, design = design)
    fit <- glmFit(dge, design = design)
    lrt <- glmLRT(fit, coef = 2)
    res_size[[paste0("Cluster_", cluster)]] <- topTags(lrt, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)$table
    
  }
  saveRDS(res_time, file = paste0(dir_out, "edgeR_time_", cell_type, "_", today, ".rds"))
  saveRDS(res_size, file = paste0(dir_out, "edgeR_size_", cell_type, "_", today, ".rds"))
}

