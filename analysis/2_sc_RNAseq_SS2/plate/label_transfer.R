rm(list = setdiff(ls(), c("so_T", "so_nT")))
today <- Sys.Date()#"2020-10-14"
dir_proj <- "/singerlab/linglin/b16/"

setwd(dir_proj)
dir_out <- paste0(dir_proj, "2_pipeline/plate/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(Matrix)
library(openxlsx)
library(ggplot2)
source(paste0(dir_proj, "1_code/make_figures/utils.R"))

## 10x data
prep_data_date <- "2020-08-06"
if (! ("so_T" %in% ls())){
  so_T <- readRDS(paste0("/singerlab/linglin/b16/2_pipeline/prep_data/so_T_", prep_data_date, ".rds"))
}
ann_T <- read.xlsx("/singerlab/linglin/b16/0_data/cluster_markers_20210826.xlsx", sheet = "T_annotation", colNames = F)
so_T$annotation <- ann_T$X2[match(paste0("C", so_T$louvain), ann_T$X1)]
so_T <- ScaleData(so_T)
so_T <- FindVariableFeatures(so_T)
so_T <- RunPCA(so_T)

## plate data
tpm_cd4 <- read.table("/singerlab/kathleen/b16/plate data/cd4/tpm_merged_cd4_20170530.selected_cells.txt") %>% as.sparse()
tpm_cd8 <- read.table("/singerlab/kathleen/b16/plate data/cd8/tpm_merged_cd8_20170526.selected_cells.txt") %>% as.sparse()
stopifnot(all(toupper(rownames(tpm_cd4))==toupper(rownames(tpm_cd8))))
rownames(tpm_cd4) <- rownames(tpm_cd8)
so_cd4 <- CreateSeuratObject(counts = tpm_cd4, min.cells = 1); rm(tpm_cd4)
so_cd8 <- CreateSeuratObject(counts = tpm_cd8, min.cells = 1); rm(tpm_cd8)
so_cd4$orig.ident <- "CD4"; so_cd8$orig.ident <- "CD8"
so_cd4 <- NormalizeData(so_cd4)
so_cd4 <- ScaleData(so_cd4)
so_cd4 <- FindVariableFeatures(so_cd4)
so_cd4 <- RunPCA(so_cd4)
ElbowPlot(so_cd4, ndims = 50)
so_cd8 <- NormalizeData(so_cd8)
so_cd8 <- ScaleData(so_cd8)
so_cd8 <- FindVariableFeatures(so_cd8)
so_cd8 <- RunPCA(so_cd8)
ElbowPlot(so_cd8, ndims = 50)
so_plate <- merge(so_cd4, so_cd8)
so_plate <- NormalizeData(so_plate)
so_plate <- ScaleData(so_plate)
so_plate <- FindVariableFeatures(so_plate, selection.method = "vst")
so_plate <- RunPCA(so_plate)
ElbowPlot(so_plate, ndims = 50)

## label transfer
# CD4
anchors <- FindTransferAnchors(reference = so_T, query = so_cd4, dims = 1:20)
predictions <- TransferData(anchorset = anchors, dims = 1:20, refdata = as.character(so_T$louvain))
so_cd4 <- AddMetaData(so_cd4, metadata = predictions)
so_cd4$predicted.id_annotation <- ann_T$X2[match(paste0("C", so_cd4$predicted.id), ann_T$X1)]
p_cd4 <- ggplot(so_cd4@meta.data) + 
  geom_bar(aes(x = predicted.id_annotation, fill = predicted.id_annotation), width = 0.7) +
  scale_fill_manual(values = get_cluster_colors(celltype = "T", cluster_name_type = "annotation")) +
  labs(x = "predicted cluster", fill = "", title ="CD4") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1))
# CD8
anchors <- FindTransferAnchors(reference = so_T, query = so_cd8, dims = 1:20)
predictions <- TransferData(anchorset = anchors, dims = 1:20, refdata = as.character(so_T$louvain))
so_cd8 <- AddMetaData(so_cd8, metadata = predictions)
so_cd8$predicted.id_annotation <- ann_T$X2[match(paste0("C", so_cd8$predicted.id), ann_T$X1)]
p_cd8 <- ggplot(so_cd8@meta.data) + 
  geom_bar(aes(x = predicted.id_annotation, fill = predicted.id_annotation), width = 0.7) +
  scale_fill_manual(values = get_cluster_colors(celltype = "T", cluster_name_type = "annotation")) +
  labs(x = "predicted cluster", fill = "", title ="CD8") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1))
# together
anchors <- FindTransferAnchors(reference = so_T, query = so_plate, dims = 1:20)
predictions <- TransferData(anchorset = anchors, dims = 1:20, refdata = as.character(so_T$louvain))
so_plate <- AddMetaData(so_plate, metadata = predictions)
so_plate$predicted.id_annotation <- ann_T$X2[match(paste0("C", so_plate$predicted.id), ann_T$X1)]
p_cd4_cd8 <- ggplot(so_plate@meta.data) + 
  geom_bar(aes(x = predicted.id_annotation, fill = orig.ident), width = 0.5, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("CD4" = "black", "CD8" = "red")) +
  labs(x = "predicted cluster", fill = "", title ="CD4 and CD8") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1))

pdf(paste0(dir_out, "label_transfer_", today, ".pdf"), width = 5, height = 5)
p_cd4
p_cd8
p_cd4_cd8
dev.off()

# to save space
so_cd4[["RNA"]]@scale.data <- matrix()
so_cd8[["RNA"]]@scale.data <- matrix()
so_plate[["RNA"]]@scale.data <- matrix()
saveRDS(so_cd4, file = paste0(dir_out, "so_cd4_", today, ".rds"))
saveRDS(so_cd8, file = paste0(dir_out, "so_cd8_", today, ".rds"))
saveRDS(so_plate, file = paste0(dir_out, "so_cd4_cd8_", today, ".rds"))

so_plate <- readRDS("/singerlab/linglin/b16/2_pipeline/plate/so_cd4_cd8_2021-08-26.rds")
library(Seurat)
FeaturePlot(so_plate, "Penk", order = TRUE)
DimPlot(so_plate, group.by = "predicted.id")
DimPlot(so_plate, group.by = "orig.ident")
VlnPlot(so_plate, "Penk", group.by = "predicted.id")
so_plate@active.ident <- factor(so_plate$predicted.id)
so_plate$Penk_positive <- factor(so_plate[["RNA"]]@counts["Penk",] > 0)
m_4 <- FindMarkers(so_plate, ident.1 = "TRUE", group.by = 'Penk_positive', subset.ident = "4")
m_7 <- FindMarkers(so_plate, ident.1 = "TRUE", group.by = 'Penk_positive', subset.ident = "7")

library(tidyverse)
m_4 <- m_4 %>% filter(p_val_adj < 0.1)
m_7 <- m_7 %>% filter(p_val_adj < 0.1)

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "T_4")
writeData(wb, "T_4", m_4, rowNames = TRUE)
addWorksheet(wb, "T_7")
writeData(wb, "T_7", m_7, rowNames = TRUE)
saveWorkbook(wb, "2_pipeline/plate/Penk_pos_vs_neg_T4_T7.xlsx", overwrite = TRUE)