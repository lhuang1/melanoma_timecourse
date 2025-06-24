res_time_t <- readRDS("/singerlab/linglin/b16/2_pipeline/TSA/DE_time_size/edgeR_time_T_2020-12-30.rds")
res_time_nt <- readRDS("/singerlab/linglin/b16/2_pipeline/TSA/DE_time_size/edgeR_time_nT_2020-12-30.rds")
res_size_t <- readRDS("/singerlab/linglin/b16/2_pipeline/TSA/DE_time_size/edgeR_size_T_2020-12-30.rds")
res_size_nt <- readRDS("/singerlab/linglin/b16/2_pipeline/TSA/DE_time_size/edgeR_size_nT_2020-12-30.rds")

thresh_lfc <- 0
thresh_fdr <- 0.05
df_time <- c(lapply(names(res_time_t), function(i){
  res_time_t[[i]] %>% rownames_to_column(var = "Gene") %>%
    # filter(FDR < thresh_fdr & abs(logFC) > thresh_lfc) %>%
    mutate(Cluster = sub("Cluster_", "T_", i))
}), lapply(names(res_time_nt), function(i){
  res_time_nt[[i]] %>% rownames_to_column(var = "Gene") %>%
    # filter(FDR < thresh_fdr & abs(logFC) > thresh_lfc) %>%
    mutate(Cluster = sub("Cluster_", "Non-T_", i))
})) %>% Reduce(rbind, .) %>% mutate(unique = !duplicated(Gene), cov = "Time")

df_size <- c(lapply(names(res_size_t), function(i){
  res_size_t[[i]] %>% rownames_to_column(var = "Gene") %>%
    # filter(FDR < thresh_fdr & abs(logFC) > thresh_lfc) %>%
    mutate(Cluster = sub("Cluster_", "T_", i))
}), lapply(names(res_size_nt), function(i){
  res_size_nt[[i]] %>% rownames_to_column(var = "Gene") %>%
    # filter(FDR < thresh_fdr & abs(logFC) > thresh_lfc) %>%
    mutate(Cluster = sub("Cluster_", "Non-T_", i))
})) %>% Reduce(rbind, .) %>% mutate(unique = !duplicated(Gene), cov = "Size")


df_time_size <- rbind(df_time, de_size)
df_time_size$cov <- factor(df_time_size$cov, levels = c("Time", "Size"))

saveRDS(df_time_size, paste0(dir_out, "edgeR_time_size_all_2020-12-30.rds"))