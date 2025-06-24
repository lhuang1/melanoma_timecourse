## dir_proj
dir_proj <- "~/project/melanoma_timecourse/"

## cluster colors
color_cluster_T <- c('T3' = '#33A02C', 'T1' = '#B15928', 'T7' = '#E31A1C', 'T4' = '#FB9A99', 'T5' = '#A6CEE3', 'T6' = '#1F78B4', 'T2' = '#CAB2D6', 'T8' = '#6A3D9A', 'T0' = '#FDBF6F', 'T9' = '#FF7F00')


color_cluster_nT <- c('Non-T0' = '#F8766D', 'Non-T2' = '#FF64B0', 'Non-T1' = '#00BFC4', 
                      'Non-T5' = '#00B4F0', 'Non-T8' = '#619CFF', 'Non-T9' = '#C77CFF', 'Non-T6' = '#F564E3', 'Non-T7' = '#00C08B',
                      'Non-T4' = '#B79F00', 'Non-T10' = '#DE8C00', 'Non-T3' = '#7CAE00', 'Non-T11' = '#00BA38', 'Non-T12' = '#999999')

### cluster annotation
wb <- loadWorkbook(paste0(dir_proj, "1_analysis/timecourse_sc_RNAseq_10x/1_clustering/cluster_markers_20241206.xlsx"))
ann_T <- readWorkbook(wb, sheet = "T_annotation", colNames = T)
ann_nT <- readWorkbook(wb, sheet = "nonT_annotation", colNames = T)
ann_nT$annotation[ann_nT$annotation == " "] <- "Other"
color_annotation_T <- color_cluster_T %>% `names<-`(ann_T$annotation[match(names(color_cluster_T), ann_T$cluster)])
color_annotation_nT <- color_cluster_nT %>% `names<-`(ann_nT$annotation[match(names(color_cluster_nT), ann_nT$cluster)])

color_cluster_annotation_T <- color_cluster_T %>% `names<-`(paste(names(color_cluster_T), names(color_annotation_T), sep = ". "))
color_cluster_annotation_nT <- color_cluster_nT %>% `names<-`(paste(names(color_cluster_nT), names(color_annotation_nT), sep = ". "))


## color vector to sample from
library(RColorBrewer)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# length(col_vector) #74



###########
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

map_signif_level_func <- function(p) {
  if (p > 0.001) return(paste0("p = ", sprintf("%.3f", p)))
  return("p < 0.001")
}

plot_feature <- function(so, g, pt_size = 0.1, min_quantile = 0, max_quantile = 1, sort_by_expr = TRUE) {
  if (g %in% colnames(so@meta.data)) {
    df <- so@meta.data[,c("UMAP_1", "UMAP_2", g)]
    color_bar_name <- "Signature\nScore"
  } else if (g %in% rownames(so)) {
    df <- cbind(so@meta.data[,c("UMAP_1", "UMAP_2")], so[["RNA"]]@data[g,])
    colnames(df)[3] <- g
    color_bar_name <- "Expression"
  } else {
    cat("\t!! Error: Feature", g, "does not exist.\n")
    return(NULL)
  }
  if (sort_by_expr) {
    df <- df[order(df[,g], decreasing = FALSE),]
  }
  
  lim_max <- quantile(df[,g], max_quantile, na.rm = T)
  lim_min <- quantile(df[,g], min_quantile, na.rm = T)
  df[,g] <- pmax(pmin(df[,g], lim_max), lim_min)
  
  p <- ggplot(df) +
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = get(g)), size = pt_size) +
    scale_color_gradient(low = "grey80", high = "blue", na.value = "grey40") +
    # scale_color_viridis_c() +
    labs(color = color_bar_name) +
    ggtitle(g) +
    theme(legend.position = "right",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(face = "bold", hjust = 0.5))
  return(p)
}




## compare marker and cluster DE results
compare_marker_deg <- function(marker_df, de_res, de_method) {
  library(ggrepel)
  cluster_vec <- sort(unique(marker_df$cluster))
  setdiff_list <- sapply(cluster_vec, function(i) {
    if (de_method == "wilcoxon") {
      up <- de_res %>% 
        filter(p_val_adj < 0.05 & avg_logFC > 0.25 & cluster == as.numeric(sub("C", "", i))) %>% 
        select(gene) %>% 
        unlist() %>% 
        as.character()
    } else if (de_method == "psb") {
      up <- de_res[[sub("C", "", i)]] %>% 
        rownames_to_column(var = "gene") %>% 
        filter(FDR < 0.05 & logFC > 0.25) %>% 
        select(gene) %>% 
        unlist() %>% 
        as.character()
    }
    setdiff(toupper(marker_df$gene[marker_df$cluster == i]), toupper(up))
  })
  return(setdiff_list)
}


plot_volcano_marker_deg <- function(marker_df, de_res, de_method, cell_type) {
  library(ggrepel)
  cluster_vec <- sort(unique(marker_df$cluster))
  p_volcano <- list()
  for (i in cluster_vec) {
    m_i <- marker_df$gene_name[marker_df$cluster == i]
    
    if (de_method == "wilcoxon") {
      tmp <- de_res[de_res$cluster == as.numeric(sub("C", "", i)) & de_res$avg_logFC > 0,]
      to_show_label <- union(m_i, tmp$gene[c(order(tmp$avg_logFC, decreasing = T)[1:20], order(tmp$p_val, decreasing = F)[1:20])])
      plt_df <- tmp %>%
        mutate(is_marker = gene %in% m_i) %>%
        mutate(label = ifelse(gene %in% to_show_label, gene, NA)) %>% 
        select(logFC = avg_logFC, FDR = p_val_adj, is_marker, label)
      x_intercepts <- c(0.25, log2(1.5))
    } else if (de_method %in% c("psb")) {
      tmp <- de_res[[sub("C", "", i)]]
      to_show_label <- union(m_i, rownames(tmp)[c(order(tmp$logFC, decreasing = T)[1:20], order(tmp$PValue, decreasing = F)[1:20])])
      plt_df <- tmp %>%
        tibble::rownames_to_column(var = "gene") %>%
        mutate(is_marker = gene %in% m_i) %>%
        mutate(label = ifelse(gene %in% to_show_label, gene, NA))
      x_intercepts <- c(-log2(1.5), -0.25, 0.25, log2(1.5))
    }
    
    p_volcano[[i]] <- plt_df %>%
      ggplot(aes(x = logFC, y = -log10(FDR), color = is_marker)) +
      geom_point() +
      geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
      geom_vline(xintercept = x_intercepts, linetype = 2, color = "grey50") +
      geom_label_repel(aes(label = label)) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
      ggtitle(paste0(cell_type, i)) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  }
  return(p_volcano)
}







plot_dot <- function(so, features, group.by, scale_by_cluster, wrap_n = 20) {
  library(stringr)
  group_vec <- levels(factor(so@meta.data[,group.by])) %>% setdiff(., "unknown")
  avg_expr <- sapply(group_vec, function(g){
    Matrix::rowMeans(so[["RNA"]]@data[features, so@meta.data[,group.by] == g])
  })
  if (scale_by_cluster) {
    avg_expr <- scale(t(avg_expr)) %>% t()
    lab_color <- "Scaled\nAverage Expression"
  } else {
    lab_color <- "Average Expression"
  }
  pct_det <- sapply(group_vec, function(g){
    Matrix::rowMeans(so[["RNA"]]@data[features, so@meta.data[,group.by] == g] > 0) * 100
  })
  plt_df <- cbind(
    gather(cbind(data.frame(avg_expr, check.names = F), gene = rownames(avg_expr)), key = "cluster", value = "avg_expr", -gene),
    gather(data.frame(pct_det, check.names = F), key = "cluster2", value = "pct_det")
  ) %>% 
    mutate(gene = factor(gene, levels = unique(features)),
           cluster = factor(str_wrap(cluster, wrap_n), levels = rev(str_wrap(group_vec, wrap_n))))
  square_trans <- scales::trans_new("square", transform = function(x){x^2}, inverse = function(x){sqrt(x)})
  p <- plt_df %>% ggplot() +
    geom_point(aes(x = gene, y = cluster, fill = avg_expr, size = pct_det), color = "black", shape = 21) +
    # scale_color_viridis_c() +
    # scale_color_gradient(low = "grey80", high = "blue") +
    scale_fill_gradient2(low = "magenta", mid = "black", high = "yellow", breaks = seq(-2.5, 2.5, length.out = 5)) +
    scale_size_continuous(range = c(0, 5), trans = square_trans, breaks = seq(0, 100, 25), limits = c(0, 100)) +
    labs(x = "", y = "", fill = lab_color, size = "Percent Expression") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)
}


get_cluster_colors <- function(celltype = "T", cluster_name_type = "name") {
  if (celltype == "nT") {
    cluster_number_nT <- unique(sort(so_nT$louvain))
    color_cluster_number_nT <- gg_color_hue(length(cluster_number_nT)); names(color_cluster_number_nT) <- as.character(cluster_number_nT)
    if (cluster_name_type == "number") return(color_cluster_number_nT)
    if (cluster_name_type == "name") {
      names(color_cluster_number_nT) <- paste0("Non-T_", names(color_cluster_number_nT))
      return(color_cluster_number_nT)
    }
    if (cluster_name_type == "annotation"){
      cluster_annotation_nT <- unique(so_nT$annotation[order(so_nT$louvain)]) %>% as.character()
      color_cluster_annotation_nT <- color_cluster_number_nT; names(color_cluster_annotation_nT) <- cluster_annotation_nT
      return(color_cluster_annotation_nT)
    }
  } else if (celltype == "T") {
    cluster_number_T <- unique(sort(so_T$louvain))
    color_cluster_number_T <- gg_color_hue(length(cluster_number_T)); names(color_cluster_number_T) <- as.character(cluster_number_T)
    if (cluster_name_type == "number") return(color_cluster_number_T)
    if (cluster_name_type == "name") {
      names(color_cluster_number_T) <- paste0("Non-T_", names(color_cluster_number_T))
      return(color_cluster_number_T)
    }
    if (cluster_name_type == "annotation"){
      cluster_annotation_T <- unique(so_T$annotation[order(so_T$louvain)]) %>% as.character()
      color_cluster_annotation_T <- color_cluster_number_T; names(color_cluster_annotation_T) <- cluster_annotation_T
      return(color_cluster_annotation_T)
    }
  }
}
