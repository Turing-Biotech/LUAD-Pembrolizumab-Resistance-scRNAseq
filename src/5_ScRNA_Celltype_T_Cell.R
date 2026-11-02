# Set ####
options(timeout = 7355608); getOption("timeout")
options(stringsAsFactors = FALSE)
getwd(); rm(list = ls()); gc()
set.seed(19981108)
if (!dir.exists("5_ScRNA_Celltype_T_Cell")) {dir.create("5_ScRNA_Celltype_T_Cell")}

# Library ####
library(Seurat); packageVersion("Seurat")
library(tidyverse); packageVersion("tidyverse")
library(harmony); packageVersion("harmony")
library(patchwork); packageVersion("patchwork")
library(ggraph); packageVersion("ggraph")
library(cowplot); packageVersion("cowplot")
library(clustree); packageVersion("clustree")
library(ggrepel); packageVersion("ggrepel")
library(ggrastr); packageVersion("ggrastr")
library(ggalluvial); packageVersion("ggalluvial")
library(ggradar2); packageVersion("ggradar2")
library(ggpubr); packageVersion("ggpubr")
library(cols4all); packageVersion("cols4all")

# Colours ####
Colours_ScRNA_Sample <- c(
  "#FF1C1CFF", "#FB8B05FF", "#FFFF00FF", "#33CC33FF", "#6699FFFF", "#1CFFD2FF", "#A246FFFF", "#FF7AC2FF")
Colours_ScRNA_Clustree <- c(
  "#FF6347FF", "#FB8B05FF", "#FFFF00FF", "#FCD337FF", "#66FF66FF", "#33CC33FF", 
  "#FF1C1CFF", "#1CCEFFFF", "#6699FFFF", "#1CFFD2FF", "#A246FFFF", "#FF7AC2FF")
Colours_ScRNA_Cell <- c(
  "#FFEB4EFF", "#9933CCFF", "#33CC33FF", "#1CCEFFFF", 
  "#FF7AC2FF", "#FF6354FF", "#66FF66FF", "#C4A1FFFF", 
  "#FFC6FFFF", "#FFA845FF")
Colours_ScRNA_Cluster <- c(
  "#FF1C1CFF", "#FFA845FF", "#FFEB4EFF", "#7CFC00FF", "#1CFFD2FF", 
  "#9BF6FFFF", "#C4A1FFFF", "#FF7AC2FF", "#DC143CFF", "#FFA845FF", 
  "#FF6354FF", "#1CCEFFFF", "#A246FFFF", "#CD36FFFF", "#CCCC66FF", 
  "#1CFF86FF", "#A0C4FFFF", "#993333FF", "#FFADADFF", "#FFD6A5FF", 
  "#FDFFB6FF", "#CAFFBFFF", "#F08080FF", "#808000FF", "#088066FF", 
  "#D2691EFF", "#FFB775FF", "#4682B4FF", "#1E90FFFF", "#6666FFFF", 
  "#5047FFFF", "#BDB2FFFF", "#FFC6FFFF", "#228B22FF", "#7B68EEFF", 
  "#D2B48CFF", "#DEB887FF", "#33CC33FF", "#33CCCCFF", "#E57373FF", 
  "#FF1493FF", "#A0522DFF", "#D595A7FF", "#6699FFFF", "#CCFFFFFF", 
  "#1B1919FF", "#6F286AFF", "#8F7700FF", "#99991EFF", "#7F8C8DFF")

# Load ####
load(file = "5_ScRNA_Celltype_Major/MetaData_Celltype_Major.RData") # MetaData_Celltype_Major
load(file = "1_ScRNA_SeuratData/Seurat_QualityControl.RData") # Seurat_QualityControl.RData
head(MetaData_Celltype_Major); head(Seurat_QualityControl@meta.data)
class(MetaData_Celltype_Major); class(Seurat_QualityControl@meta.data)
dim(MetaData_Celltype_Major); dim(Seurat_QualityControl@meta.data)
Seurat_QualityControl@meta.data <- Seurat_QualityControl@meta.data %>%
  tibble::rownames_to_column("Cell_ID") %>%
  dplyr::left_join(MetaData_Celltype_Major %>%
                     tibble::rownames_to_column("Cell_ID") %>%
                     dplyr::select(Cell_ID, Resolution_1, Celltype_Major), by = "Cell_ID") %>%
  dplyr::rename(Resolution_Celltype_Major = Resolution_1) %>% 
  tibble::column_to_rownames("Cell_ID") %>% 
  base::as.data.frame(); head(Seurat_QualityControl@meta.data); table(Seurat_QualityControl$Celltype_Major)

# Subpopulation Extraction ####
Seurat_Celltype_T_Cell <- subset(Seurat_QualityControl, subset = Celltype_Major == "T_Cell")
head(Seurat_Celltype_T_Cell@meta.data); Seurat_Celltype_T_Cell; rm(Seurat_QualityControl, MetaData_Celltype_Major); gc()

# Integration NormalizeData => FindVariableFeatures => ScaleData ####
# NormalizeData
Seurat_Celltype_T_Cell <- JoinLayers(object = Seurat_Celltype_T_Cell)
Seurat_Celltype_T_Cell <- NormalizeData(object = Seurat_Celltype_T_Cell,
                                           normalization.method = "LogNormalize", 
                                           scale.factor = 10000,
                                           verbose = TRUE)
# FindVariableFeatures
Seurat_Celltype_T_Cell <- FindVariableFeatures(object = Seurat_Celltype_T_Cell, 
                                                  selection.method = "vst", 
                                                  nfeatures = 3000, # [!]
                                                  verbose = TRUE)

# CellCycle
S_Genes_Human <- Seurat::cc.genes.updated.2019$s.genes
G2M_Genes_Human <- Seurat::cc.genes.updated.2019$g2m.genes
S_Genes_Match <- CaseMatch(S_Genes_Human, rownames(Seurat_Celltype_T_Cell))
G2M_Genes_Match <- CaseMatch(G2M_Genes_Human, rownames(Seurat_Celltype_T_Cell))
Seurat_Celltype_T_Cell <- CellCycleScoring(object = Seurat_Celltype_T_Cell,
                                           s.features = S_Genes_Match,
                                           g2m.features = G2M_Genes_Match,
                                           set.ident = FALSE)
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
table(Seurat_Celltype_T_Cell@meta.data$Phase)

# ScaleData
Seurat_Celltype_T_Cell <- ScaleData(object = Seurat_Celltype_T_Cell, 
                                       vars.to.regress = c("S.Score", "G2M.Score", "Percent_MT"),
                                       features = VariableFeatures(object = Seurat_Celltype_T_Cell)); Seurat_Celltype_T_Cell
DefaultAssay(Seurat_Celltype_T_Cell)

# SCT ####
Seurat_Celltype_T_Cell <- SCTransform(object = Seurat_Celltype_T_Cell, 
                                         assay = "RNA",
                                         new.assay.name = "SCT",
                                         variable.features.n = 3000, # [!]             
                                         vars.to.regress = c("S.Score", "G2M.Score", "Percent_MT"),
                                         verbose = TRUE)
DefaultAssay(Seurat_Celltype_T_Cell)

# PCA ####
Seurat_Celltype_T_Cell <- RunPCA(object = Seurat_Celltype_T_Cell, 
                                    reduction.name = "pca",
                                    features = VariableFeatures(object = Seurat_Celltype_T_Cell),
                                    npcs = 50, 
                                    seed.use = 19981108,
                                    verbose = TRUE)
names(Seurat_Celltype_T_Cell@reductions)
ElbowPlot(Seurat_Celltype_T_Cell, ndims = 50, reduction = "pca")
Seurat_Celltype_T_Cell

# Dims
PC_Stdev <- Seurat_Celltype_T_Cell@reductions$pca@stdev
PC_Explanation_Rate <- PC_Stdev^2 / sum(PC_Stdev^2) * 100
PC_Explanation_Rate
PC_Explanation_Rate_Cumulative <- cumsum(PC_Explanation_Rate)
PC_Explanation_Rate_Cumulative[1:50] 
PC_Explanation_Rate_Cumulative[1:40] 
PC_Explanation_Rate_Cumulative[1:30] 
PC_Explanation_Rate_Cumulative[1:20] 
PC_Explanation_Rate_Cumulative[1:10]
Dims <- 1:30

# FindNeighbors & FindClusters
names(Seurat_Celltype_T_Cell@reductions)
Seurat_Celltype_T_Cell <- FindNeighbors(object = Seurat_Celltype_T_Cell, 
                                        reduction = "pca", 
                                        dims = Dims)
Seurat_Celltype_T_Cell <- FindClusters(object = Seurat_Celltype_T_Cell, 
                                       resolution = 1, 
                                       cluster.name = "PCA_1",
                                       random.seed = 19981108,
                                       algorithm = 4)
head(Seurat_Celltype_T_Cell@meta.data, n = 10)

# UMAP & TSNE [PCA]
# UMAP
Seurat_Celltype_T_Cell <- RunUMAP(object = Seurat_Celltype_T_Cell, 
                                     dims = Dims,
                                     min.dist = 0.3, # 0.3 [-]
                                     n.neighbors = 50L, # 30L [+]
                                     spread = 0.6, # 1 [-]
                                     reduction.name = "pca_umap",
                                     reduction = "pca")
# TSNE
Seurat_Celltype_T_Cell <- RunTSNE(object = Seurat_Celltype_T_Cell,
                                     dims = Dims,
                                     seed.use = 19981108,
                                     perplexity = 30,
                                     reduction = "pca",
                                     reduction.name = "pca_tsne",
                                     verbose = TRUE)
names(Seurat_Celltype_T_Cell@reductions)

# Plot_PCA_Sample
Plot_PCA_Sample <- DimPlot(object = Seurat_Celltype_T_Cell, 
                           reduction = "pca", 
                           group.by = "Sample", 
                           pt.size = 0.5, 
                           seed = 19981108,
                           raster = TRUE) +
  labs(title = "PCA", 
       x = "PC 1", 
       y = "PC 2") +  
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_manual(values = Colours_ScRNA_Sample)
Plot_PCA_Sample

# Plot_UMAP:PCA
Plot_UMAP_PCA <- DimPlot(object = Seurat_Celltype_T_Cell, 
                         reduction = "pca_umap",
                         label = FALSE,
                         group.by = "Sample", 
                         raster = TRUE, 
                         pt.size = 0.5,
                         label.size = 5) + 
  theme_bw()+
  labs(title = "UMAP (PCA)", 
       x = "UMAP_PCA 1", 
       y = "UMAP_PCA 2") +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_manual(values = Colours_ScRNA_Sample)
Plot_UMAP_PCA

# Plot_TSNE:PCA
Plot_TSNE_PCA <- DimPlot(object = Seurat_Celltype_T_Cell, 
                         reduction = "pca_tsne",
                         label = FALSE,
                         group.by = "Sample", 
                         raster = TRUE, 
                         pt.size = 0.5,
                         label.size = 5) + 
  theme_bw()+
  labs(title = "TSNE (PCA)", 
       x = "TSNE_PCA 1", 
       y = "TSNE_PCA 2") +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_manual(values = Colours_ScRNA_Sample)
Plot_TSNE_PCA

# Harmony ####
names(Seurat_Celltype_T_Cell@reductions)
Seurat_Celltype_T_Cell <- RunHarmony(Seurat_Celltype_T_Cell,
                                        group.by.vars = "orig.ident",
                                        reduction.use = "pca",
                                        dims = Dims, 
                                        reduction.save = "harmony",
                                        project.dim = FALSE,
                                        plot_convergence = TRUE)

# FindNeighbors & FindClusters
Seurat_Celltype_T_Cell <- FindNeighbors(object = Seurat_Celltype_T_Cell, 
                                           reduction = "harmony", 
                                           dims = Dims)
Seurat_Celltype_T_Cell <- FindClusters(object = Seurat_Celltype_T_Cell, 
                                          resolution = 1, 
                                          cluster.name = "Harmony_1",
                                          random.seed = 19981108,
                                          algorithm = 4)
print(names(Seurat_Celltype_T_Cell@reductions))

# UMAP & TSNE [Harmony]
Seurat_Celltype_T_Cell <- RunUMAP(object = Seurat_Celltype_T_Cell, 
                                     dims = Dims,
                                     min.dist = 0.3, # 0.3 [-]
                                     n.neighbors = 50L, # 30L [+]
                                     spread = 0.6, # 1 [-]
                                     reduction.name = "harmony_umap",
                                     reduction = "harmony")
Seurat_Celltype_T_Cell <- RunTSNE(object = Seurat_Celltype_T_Cell,
                                     dims = Dims,
                                     seed.use = 19981108,
                                     perplexity = 30,
                                     reduction = "harmony",
                                     reduction.name = "harmony_tsne",
                                     verbose = TRUE)
print(names(Seurat_Celltype_T_Cell@reductions))

# Plot_Harmony_Sample
Plot_Harmony_Sample <- DimPlot(object = Seurat_Celltype_T_Cell, 
                               reduction = "harmony", 
                               group.by = "Sample", 
                               pt.size = 0.5, 
                               seed = 19981108,
                               raster = TRUE) +
  labs(title = "Harmony", 
       x = "Harmony 1", 
       y = "Harmony 2") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_manual(values = Colours_ScRNA_Sample)
Plot_Harmony_Sample

# Plot_UMAP:Harmony
Plot_UMAP_Harmony <- DimPlot(object = Seurat_Celltype_T_Cell, 
                             reduction = "harmony_umap",
                             label = FALSE,
                             group.by = "Sample", 
                             raster = TRUE, 
                             pt.size = 0.5,
                             label.size = 5) + 
  theme_bw()+
  labs(title = "UMAP (Harmony)", 
       x = "UMAP_Harmony 1", 
       y = "UMAP_Harmony 2") +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_manual(values = Colours_ScRNA_Sample)
Plot_UMAP_Harmony

# Plot_TSNE:Harmony
Plot_TSNE_Harmony <- DimPlot(object = Seurat_Celltype_T_Cell, 
                             reduction = "harmony_tsne",
                             label = FALSE,
                             group.by = "Sample", 
                             raster = TRUE, 
                             pt.size = 0.5,
                             label.size = 5) + 
  theme_bw()+
  labs(title = "TSNE (Harmony)", 
       x = "TSNE_Harmony 1", 
       y = "TSNE_Harmony 2") +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_manual(values = Colours_ScRNA_Sample)
Plot_TSNE_Harmony

# Patchwork ####
Plot_UMAP_TSNE_PCA <- Plot_PCA_Sample | Plot_UMAP_PCA | Plot_TSNE_PCA
Plot_UMAP_TSNE_PCA
Plot_UMAP_TSNE_Harmony <- Plot_Harmony_Sample | Plot_UMAP_Harmony | Plot_TSNE_Harmony
Plot_UMAP_TSNE_Harmony
Plot_All_Integration <- 
  (Plot_PCA_Sample | Plot_UMAP_PCA | Plot_TSNE_PCA) /
  (Plot_Harmony_Sample | Plot_UMAP_Harmony | Plot_TSNE_Harmony)
Plot_All_Integration
ggsave(plot = Plot_All_Integration, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_All_Integration.pdf",
       width = 15, height = 8, units = "in", dpi = 300)

# FindNeighbors & FindClusters ####
names(Seurat_Celltype_T_Cell@reductions)
# Seurat_Celltype_T_Cell <- FindNeighbors(object = Seurat_Celltype_T_Cell, reduction = "harmony", dims = Dims, k.param = 20) # Harmony
Seurat_Celltype_T_Cell <- FindNeighbors(object = Seurat_Celltype_T_Cell, reduction = "harmony", dims = Dims, k.param = 10) # Harmony

# FindClusters
Resolution <- c(0.01, 0.05, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)
for (Res in Resolution) {
  Seurat_Celltype_T_Cell <- FindClusters(object = Seurat_Celltype_T_Cell, 
                                            resolution = Res, 
                                            cluster.name = paste0("Resolution_", Res),
                                            random.seed = 19981108,
                                            algorithm = 4)}
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
for (Res in Resolution) {
  colname <- paste0("Resolution_", Res)
  cat("Resolution:", Res, "\n")
  print(table(Seurat_Celltype_T_Cell[[colname]]))}

# UMAP & TSNE
Seurat_Celltype_T_Cell <- RunUMAP(object = Seurat_Celltype_T_Cell, 
                                     dims = Dims, 
                                     min.dist = 0.3, # 0.3 [-]
                                     n.neighbors = 50L, # 30L [+]
                                     spread = 0.6, # 1 [-]
                                     reduction.name = "harmony_umap_cluster", # !!!
                                     reduction = "harmony") # !!!
Seurat_Celltype_T_Cell <- RunTSNE(object = Seurat_Celltype_T_Cell,
                                     dims = Dims, 
                                     seed.use = 19981108,
                                     perplexity = 30,
                                     reduction = "harmony", # !!!
                                     reduction.name = "harmony_tsne_cluster", # !!!
                                     verbose = TRUE)
names(Seurat_Celltype_T_Cell@reductions)

# Plot UMAP & TSNE Resolution ####
names(Seurat_Celltype_T_Cell@reductions)
# Plot_UMAP_Resolution
Plot_UMAP_Resolution <- cowplot::plot_grid(
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_0.01") + ggtitle("Harmony: 0.01"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_0.05") + ggtitle("Harmony: 0.05"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_0.2") + ggtitle("Harmony: 0.2"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_0.4") + ggtitle("Harmony: 0.4"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_0.6") + ggtitle("Harmony: 0.6"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_0.8") + ggtitle("Harmony: 0.8"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_1")   + ggtitle("Harmony: 1.0"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_1.2") + ggtitle("Harmony: 1.2"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_1.4") + ggtitle("Harmony: 1.4"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_1.6") + ggtitle("Harmony: 1.6"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_1.8") + ggtitle("Harmony: 1.8"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_2") + ggtitle("Harmony: 2.0"),
  ncol = 4)
ggsave(plot = Plot_UMAP_Resolution, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_UMAP_Resolution.pdf",
       width = 20, height = 12, units = "in", dpi = 300)

# Plot_TSNE_Resolution
Plot_TSNE_Resolution <- cowplot::plot_grid(
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_tsne_cluster", label = TRUE, group.by = "Resolution_0.01") + ggtitle("Harmony: 0.01"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_tsne_cluster", label = TRUE, group.by = "Resolution_0.05") + ggtitle("Harmony: 0.05"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_tsne_cluster", label = TRUE, group.by = "Resolution_0.2") + ggtitle("Harmony: 0.2"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_tsne_cluster", label = TRUE, group.by = "Resolution_0.4") + ggtitle("Harmony: 0.4"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_tsne_cluster", label = TRUE, group.by = "Resolution_0.6") + ggtitle("Harmony: 0.6"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_tsne_cluster", label = TRUE, group.by = "Resolution_0.8") + ggtitle("Harmony: 0.8"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_tsne_cluster", label = TRUE, group.by = "Resolution_1")   + ggtitle("Harmony: 1.0"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_tsne_cluster", label = TRUE, group.by = "Resolution_1.2") + ggtitle("Harmony: 1.2"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_tsne_cluster", label = TRUE, group.by = "Resolution_1.4") + ggtitle("Harmony: 1.4"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_tsne_cluster", label = TRUE, group.by = "Resolution_1.6") + ggtitle("Harmony: 1.6"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_tsne_cluster", label = TRUE, group.by = "Resolution_1.8") + ggtitle("Harmony: 1.8"),
  DimPlot(object = Seurat_Celltype_T_Cell, raster = TRUE, reduction = "harmony_tsne_cluster", label = TRUE, group.by = "Resolution_2") + ggtitle("Harmony: 2.0"),
  ncol = 4)
ggsave(plot = Plot_TSNE_Resolution,
       filename = "5_ScRNA_Celltype_T_Cell/Plot_TSNE_Resolution.pdf",
       width = 20, height = 12, units = "in", dpi = 300)

# Plot_Clustree
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
Prefix <- Seurat_Celltype_T_Cell@meta.data %>% 
  tibble::as_tibble(rownames = "Cell_ID") %>% 
  dplyr::select(matches("^Resolution_[0-9.]+$")) %>% 
  names()
Prefix
Plot_Clustree <- Seurat_Celltype_T_Cell@meta.data %>%
  dplyr::select(all_of(Prefix)) %>%
  clustree::clustree(prefix = "Resolution_") + 
  labs(colour = "Harmony") + 
  guides(color = guide_legend(override.aes = list(size = 8))) +
  scale_color_manual(values = Colours_ScRNA_Clustree)
Plot_Clustree
ggsave(plot = Plot_Clustree, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_Clustree.pdf", 
       width = 15, height = 10, units = "in", dpi = 300)

# Resolution
Seurat_Celltype_T_Cell@meta.data <- Seurat_Celltype_T_Cell@meta.data %>%
  dplyr::mutate(dplyr::across(all_of(grep("^Resolution_", colnames(Seurat_Celltype_T_Cell@meta.data), value = TRUE)),
                              ~ factor(as.integer(as.character(.)), levels = sort(unique(as.integer(as.character(.)))))))
str(Seurat_Celltype_T_Cell@meta.data)
for (Res in Resolution) {
  colname <- paste0("Resolution_", Res)
  cat("Resolution:", Res, "\n")
  print(table(Seurat_Celltype_T_Cell[[colname]]))}

# UMAP & TSNE Coord
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
class(Seurat_Celltype_T_Cell@meta.data)
head(Seurat_Celltype_T_Cell@reductions$harmony_umap_cluster@cell.embeddings, n = 10) # harmony_umap_cluster
head(Seurat_Celltype_T_Cell@reductions$harmony_tsne_cluster@cell.embeddings, n = 10) # harmony_tsne_cluster
class(Seurat_Celltype_T_Cell@reductions$harmony_umap_cluster@cell.embeddings)
class(Seurat_Celltype_T_Cell@reductions$harmony_tsne_cluster@cell.embeddings)

# Extract
UMAP_Coord <- Seurat_Celltype_T_Cell@reductions$harmony_umap_cluster@cell.embeddings %>% 
  base::as.data.frame() %>% 
  dplyr::rename(UMAP_1 = harmonyumapcluster_1) %>% 
  dplyr::rename(UMAP_2 = harmonyumapcluster_2)
all(rownames(Seurat_Celltype_T_Cell@meta.data)  == rownames(UMAP_Coord))
which(rownames(Seurat_Celltype_T_Cell@meta.data) != rownames(UMAP_Coord))
TSNE_Coord <- Seurat_Celltype_T_Cell@reductions$harmony_tsne_cluster@cell.embeddings %>% 
  base::as.data.frame() %>% 
  dplyr::rename(TSNE_1 = harmonytsnecluster_1) %>% 
  dplyr::rename(TSNE_2 = harmonytsnecluster_2)
all(rownames(Seurat_Celltype_T_Cell@meta.data)  == rownames(TSNE_Coord))
which(rownames(Seurat_Celltype_T_Cell@meta.data) != rownames(TSNE_Coord))

# Bind
Seurat_Celltype_T_Cell@meta.data <- Seurat_Celltype_T_Cell@meta.data %>% 
  dplyr::bind_cols(UMAP_Coord) 
Seurat_Celltype_T_Cell@meta.data <- Seurat_Celltype_T_Cell@meta.data %>% 
  dplyr::bind_cols(TSNE_Coord)
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
class(Seurat_Celltype_T_Cell@meta.data)

# Save & Load ####
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
save(Seurat_Celltype_T_Cell, file = "1_ScRNA_SeuratData/Seurat_Celltype_T_Cell.RData")
# load(file = "1_ScRNA_SeuratData/Seurat_Celltype_T_Cell.RData")
Seurat_Celltype_T_Cell
DefaultAssay(Seurat_Celltype_T_Cell)
DefaultAssay(Seurat_Celltype_T_Cell) <- "RNA"
names(Seurat_Celltype_T_Cell@reductions)
names(Seurat_Celltype_T_Cell@assays$RNA@layers)
table(Seurat_Celltype_T_Cell@active.ident) 
head(Seurat_Celltype_T_Cell@meta.data, n = 10) 
Resolution_Select <- "Resolution_1" # !!!
Idents(Seurat_Celltype_T_Cell) <- Resolution_Select
table(Seurat_Celltype_T_Cell@active.ident) 
names(Seurat_Celltype_T_Cell@reductions)

# UMAP & TSNE Cluster ####
source(file = "0_Axes.R")
MetaData_Cluster <- Seurat_Celltype_T_Cell@meta.data
head(MetaData_Cluster, n = 10)
Plot_UMAP_Cluster <- ggplot(data = MetaData_Cluster %>% # !!!
                              dplyr::group_by(Resolution_1) %>% # !!!
                              dplyr::mutate(Label = paste0(Resolution_1, " (", n(), ")")) %>% # !!!
                              dplyr::mutate(Label = as.factor(Label)) %>% 
                              dplyr::ungroup(), 
                            mapping = aes(x = UMAP_1, # !!!
                                          y = UMAP_2, # !!!
                                          color = Label)) +
  geom_point_rast(size = 0, alpha = 0.9, raster.dpi = 300, scale = 0.5) + 
  geom_text_repel(data = MetaData_Cluster %>%
                    dplyr::group_by(Label = as.factor(Resolution_1)) %>% # !!!
                    dplyr::summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), # !!!
                  mapping = aes(label = Label),
                  size = 4,
                  fontface = "bold",
                  bg.color = "#FFFFFF", 
                  color = "#2E2E2E", 
                  segment.color = "#2E2E2E", 
                  segment.size = 1, 
                  max.overlaps = 1000, 
                  force = 2) +
  labs(title = "UMAP (Harmony:Resolution_1)", x = "", y = "", color = "Cluster") + # !!!
  scale_color_manual(values = Colours_ScRNA_Cluster) + # !!!
  guides(color = guide_legend(override.aes = list(size = 8))) +
  NoAxes() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
Plot_UMAP_Cluster <- Add_Left_Axes(Plot = Plot_UMAP_Cluster, Meta_Data = MetaData_Cluster, # !!! 
                                   X_Col = "UMAP_1", Y_Col = "UMAP_2", # !!!
                                   X_Label = "UMAP 1", Y_Label = "UMAP 2", # !!!
                                   X_Label_Offset = 1, Y_Label_Offset = 1); Plot_UMAP_Cluster
ggsave(plot = Plot_UMAP_Cluster, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_UMAP_Cluster.pdf", 
       width = 8, height = 6, dpi = 300, units = "in")

# TSNE Cluster
Plot_TSNE_Cluster <- ggplot(data = MetaData_Cluster %>% # !!!
                              dplyr::group_by(Resolution_1) %>% # !!!
                              dplyr::mutate(Label = paste0(Resolution_1, " (", n(), ")")) %>% # !!!
                              dplyr::mutate(Label = as.factor(Label)) %>%  # !!!
                              dplyr::ungroup(), 
                            mapping = aes(x = TSNE_1, # !!!
                                          y = TSNE_2, # !!!
                                          color = Label)) +
  geom_point_rast(size = 0, alpha = 0.9, raster.dpi = 300, scale = 0.5) + 
  geom_text_repel(data = MetaData_Cluster %>% 
                    dplyr::group_by(Label = as.factor(Resolution_1)) %>% # !!!
                    dplyr::summarize(TSNE_1 = median(TSNE_1), TSNE_2 = median(TSNE_2)), # !!!
                  mapping = aes(label = Label),
                  size = 4,
                  fontface = "bold",
                  bg.color = "#FFFFFF", 
                  color = "#2E2E2E", 
                  segment.color = "#2E2E2E", 
                  segment.size = 1, 
                  max.overlaps = 1000, 
                  force = 2) +
  labs(title = "TSNE (Harmony:Resolution_1)", x = "", y = "", color = "Cluster") + # !!!
  scale_color_manual(values = Colours_ScRNA_Cluster) + # !!!
  guides(color = guide_legend(override.aes = list(size = 8))) +
  NoAxes() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
Plot_TSNE_Cluster <- Add_Left_Axes(Plot = Plot_TSNE_Cluster, Meta_Data = MetaData_Cluster, 
                                   X_Col = "TSNE_1", Y_Col = "TSNE_2", 
                                   X_Label = "TSNE 1", Y_Label = "TSNE 2", 
                                   X_Label_Offset = 2, Y_Label_Offset = 2); Plot_TSNE_Cluster
ggsave(plot = Plot_TSNE_Cluster, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_TSNE_Cluster.pdf", 
       width = 8, height = 6, dpi = 300, units = "in")

# FindAllMarkers [Cluster] ####
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
FindAllMarkers_Cluster <- FindAllMarkers(object = Seurat_Celltype_T_Cell, 
                                         only.pos = TRUE,
                                         group.by = Resolution_Select,
                                         min.pct = 0.2,
                                         test.use = "wilcox", 
                                         logfc.threshold = 0.5,
                                         verbose = TRUE)
write_csv(x = FindAllMarkers_Cluster, file = "5_ScRNA_Celltype_T_Cell/FindAllMarkers_Cluster.csv")

# Marker Cheak ####
Markers_Cheak <- tibble::tribble(
  ~Celltype, ~Markers,
  "Naive_Marker", "CCR7",
  "Naive_Marker", "IL7R",
  "Naive_Marker", "LEF1",
  "Naive_Marker", "TCF7",
  "Naive_Marker", "CD27",
  "Naive_Marker", "CD28",

  "Memory_Marker", "FAS",
  "Memory_Marker", "ITGAL",
  "Memory_Marker", "CD58",
  "Memory_Marker", "IL2RB",
  "Memory_Marker", "PTPRC",

  "EOMES_Marker", "EOMES",

  "Effector_Marker", "GZMB",
  "Effector_Marker", "GZMA",
  "Effector_Marker", "GNLY",
  "Effector_Marker", "NKG7",
  "Effector_Marker", "PRF1",

  "Regulatory_Marker", "FOXP3",
  "Regulatory_Marker", "CTLA4",

  "CD4+_Marker", "CD4", 
  "CD8+_Marker", "CD8A", 
  "CD8+_Marker", "CD8B", 
)

# Markers Cheak [Plot_Dot_Markers_Cheak] ####
Plot_Dot_Markers_Cheak <- ggplot(data = DotPlot(object = Seurat_Celltype_T_Cell, 
                                                group.by = Resolution_Select, 
                                                # split.by = "Group", 
                                                # cols = c("#FFD6A5", "#CAFFBF", "#9BF6FF"),
                                                features = Markers_Cheak$Markers)$data %>%
                                   dplyr::rename(Percent_Expressed = pct.exp, 
                                                 Average_Expression = avg.exp.scaled) %>% 
                                   dplyr::left_join(Markers_Cheak, 
                                                    by = c("features.plot" = "Markers")), 
                                 mapping = aes(x = features.plot, 
                                               y = id, 
                                               size = Percent_Expressed, 
                                               fill = Average_Expression)) +
  geom_point(shape = 21, color = "#000000FF") +
  scale_size_continuous(range = c(0, 10)) +
  facet_grid(. ~ Celltype, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#5DB1DDFF", high = "#CC0000FF", midpoint = 0) +
  labs(title = "Marker Gene Expression in Cluster", x = "", y = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        strip.text = element_text(size = 10, face = "bold", angle = 90, hjust = 0), # angle = 90, hjust = 0 => 垂直
        strip.background = element_rect(fill = "#FFFFFFFF", color = "#000000FF"),
        legend.position = "right"); Plot_Dot_Markers_Cheak
ggsave(plot = Plot_Dot_Markers_Cheak, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_Dot_Markers_Cheak.pdf", 
       width = 13, height = 10, dpi = 300, units = "in")

# Marker Annotation ####
# Resolution_1 %in% c("", "") ~ "", 
Seurat_Celltype_T_Cell@meta.data <- Seurat_Celltype_T_Cell@meta.data %>%
  dplyr::mutate(Celltype_T_Cell = dplyr::case_when(
    Resolution_1 %in% c("1", "5", "8", "9", "10", "20", "21", "15") ~ "CD8+_Tem", 
    Resolution_1 %in% c("7") ~ "CD4+_Tem", 
    Resolution_1 %in% c("2", "4", "14", "16", "17", "19") ~ "CD4+_Tna", 
    Resolution_1 %in% c("18") ~ "CD8+_Tna", 
    Resolution_1 %in% c("3") ~ "CD4+_Tscm", 
    Resolution_1 %in% c("6", "11", "13") ~ "CD4+_Tcm", 
    Resolution_1 %in% c("12") ~ "CD4+_Treg", 
    TRUE ~ as.character(Resolution_1)))
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
table(Seurat_Celltype_T_Cell@meta.data$Celltype_T_Cell)

# Test ####
Plot_UMAP_TEST <- DimPlot(object = Seurat_Celltype_T_Cell, 
                          reduction = "harmony_umap_cluster", # !!!
                          group.by = "Celltype_T_Cell",  
                          label = TRUE, 
                          label.size = 5,  
                          repel = TRUE, 
                          raster = TRUE,
                          pt.size = 0.5) +  
  labs(title = "UMAP by Celltype [Test]", 
       x = "UMAP 1", 
       y = "UMAP 2") +  
  theme_bw() +
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
  scale_color_manual(values = Colours_ScRNA_Cluster) +
  scale_x_continuous(breaks = seq(-20, 20, by = 1)) +
  scale_y_continuous(breaks = seq(-20, 20, by = 1)); Plot_UMAP_TEST
ggsave(plot = Plot_UMAP_TEST, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_UMAP_TEST.pdf", 
       width = 8, height = 6, dpi = 300, units = "in")

# Marker Celltype ####
Markers_Celltype <- tibble::tribble(
  ~Celltype, ~Markers,
  "Naive_Marker", "CCR7",
  "Naive_Marker", "IL7R",
  "Naive_Marker", "LEF1",
  "Naive_Marker", "TCF7",
  "Naive_Marker", "CD27",
  "Naive_Marker", "CD28",

  "Memory_Marker", "FAS",
  "Memory_Marker", "ITGAL",
  "Memory_Marker", "CD58",
  "Memory_Marker", "IL2RB",
  "Memory_Marker", "PTPRC",

  "EOMES_Marker", "EOMES",

  "Effector_Marker", "GZMB",
  "Effector_Marker", "GZMA",
  "Effector_Marker", "GNLY",
  "Effector_Marker", "NKG7",
  "Effector_Marker", "PRF1",

  "Regulatory_Marker", "FOXP3",
  "Regulatory_Marker", "CTLA4",

  "CD4+_Marker", "CD4", 
  "CD8+_Marker", "CD8B", 
  )

# Markers Celltype [Plot_Dot_Markers_Celltype] ####
Plot_Dot_Markers_Celltype <- ggplot(data = DotPlot(object = Seurat_Celltype_T_Cell, 
                                                   group.by = "Celltype_T_Cell", 
                                                   # split.by = "Group", 
                                                   features = Markers_Celltype$Markers)$data %>%
                                      dplyr::left_join(Markers_Celltype, 
                                                       by = c("features.plot" = "Markers")), 
                                    mapping = aes(x = features.plot, 
                                                  y = id, 
                                                  size = pct.exp, 
                                                  fill = avg.exp.scaled)) +
  geom_point(shape = 21, color = "#000000FF") +
  scale_size_continuous(range = c(0, 8), name = "% Expressed", breaks = c(25, 50, 75, 100), limits = c(0, 100)) +
  facet_grid(. ~ Celltype, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#5DB1DDFF", high = "#CC0000FF", midpoint = 0, name = "Expression") +
  theme_classic() +
  labs(title = "Marker Gene Expression in Celltype", x = "", y = "") +
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.grid = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold", angle = 90, hjust = 0), # angle = 90, hjust = 0 => 垂直
        strip.background = element_rect(fill = "#FFFFFFFF", color = "#FFFFFFFF"), 
        legend.position = "right"); Plot_Dot_Markers_Celltype
ggsave(plot = Plot_Dot_Markers_Celltype, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_Dot_Markers_Celltype.pdf", 
       width = 12, height = 5, dpi = 300, units = "in")

# FindAllMarkers [Celltype] ####
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
FindAllMarkers_Celltype <- FindAllMarkers(object = Seurat_Celltype_T_Cell, 
                                          only.pos = TRUE,
                                          group.by = "Celltype_T_Cell",
                                          min.pct = 0.2,
                                          test.use = "wilcox", 
                                          logfc.threshold = 0.5,
                                          verbose = TRUE)
write_csv(x = FindAllMarkers_Celltype, file = "5_ScRNA_Celltype_T_Cell/FindAllMarkers_Celltype.csv")

# UMAP & TSNE Celltype ####
source(file = "0_Axes.R")
MetaData_Celltype_T_Cell <- Seurat_Celltype_T_Cell@meta.data
MetaData_Celltype_T_Cell
save(MetaData_Celltype_T_Cell, file = "5_ScRNA_Celltype_T_Cell/MetaData_Celltype_T_Cell.RData")

# UMAP Celltype
Plot_UMAP_Celltype <- ggplot(data = MetaData_Celltype_T_Cell %>% # !!!
                               dplyr::group_by(Celltype_T_Cell) %>% # !!!
                               dplyr::mutate(Label = paste0(Celltype_T_Cell, " (", n(), ")")) %>% # !!!
                               dplyr::mutate(Label = as.factor(Label)) %>% 
                               dplyr::ungroup(), 
                             mapping = aes(x = UMAP_1, # !!!
                                           y = UMAP_2, # !!!
                                           color = Label)) +
  geom_point_rast(size = 0, alpha = 0.8, raster.dpi = 300, scale = 0.5) + 
  geom_text_repel(data = MetaData_Celltype_T_Cell %>%
                    dplyr::group_by(Label = as.factor(Celltype_T_Cell)) %>% # !!!
                    dplyr::summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), # !!!
                  mapping = aes(label = Label),
                  size = 4,
                  fontface = "bold",
                  bg.color = "#FFFFFFFF", 
                  color = "#2E2E2EFF", 
                  segment.color = "#2E2E2EFF", 
                  segment.size = 1, 
                  max.overlaps = 1000, 
                  force = 2) +
  labs(title = "UMAP (Celltype_T_Cell)", x = "", y = "", color = "Celltype_T_Cell") + # !!!
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = Colours_ScRNA_Cell) + # !!!
  NoAxes() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
Plot_UMAP_Celltype <- Add_Left_Axes(Plot = Plot_UMAP_Celltype, Meta_Data = MetaData_Celltype_T_Cell, 
                                    X_Col = "UMAP_1", Y_Col = "UMAP_2", 
                                    X_Label = "UMAP 1", Y_Label = "UMAP 2", 
                                    X_Label_Offset = 1, Y_Label_Offset = 1); Plot_UMAP_Celltype
ggsave(plot = Plot_UMAP_Celltype, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_UMAP_Celltype.pdf", 
       width = 8, height = 6, dpi = 300, units = "in")

# Save ####
save(Seurat_Celltype_T_Cell, file = "1_ScRNA_SeuratData/Seurat_Celltype_T_Cell.RData")

# Cell Ratio Stacked Bar Plot ####
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
Plot_Celltype_Sample <- Seurat_Celltype_T_Cell@meta.data %>%
  dplyr::group_by(Sample, Celltype_T_Cell) %>%
  dplyr::summarise(Count = n(), .groups = "drop_last") %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup() %>% 
  ggplot2::ggplot(mapping = aes(x = Sample, y = Proportion, 
                                fill = Celltype_T_Cell,
                                stratum = Celltype_T_Cell, 
                                alluvium = Celltype_T_Cell)) +
  ggplot2::scale_y_continuous(labels = scales::percent) +  
  ggplot2::scale_fill_manual(values = Colours_ScRNA_Cell) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Proportion of Cell Types (Sample)",
                x = "Sample", y = "Proportion", fill = "Cell Type") +
  ggplot2::theme(panel.grid = element_blank(), 
                 axis.title.x = element_text(size = 12, face = "bold"), 
                 axis.title.y = element_text(size = 12, face = "bold"), 
                 axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
                 axis.text.y = element_text(size = 10, face = "bold"), 
                 legend.text = element_text(size = 10, face = "bold"), 
                 legend.title = element_text(size = 12, face = "bold"), 
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                 legend.key.height = unit(3, "lines"),
                 legend.key.width  = unit(3, "lines")) + 
  ggplot2::geom_col(width = 0.5, color = NA, linewidth = 1) + 
  ggalluvial::geom_flow(width = 0.5, alpha = 0.25, knot.pos = 0.25, color = "#FFFFFFFF", linewidth = 0.5) + 
  ggalluvial::geom_alluvium(width = 0.5, alpha = 1, knot.pos = 0.25, fill = NA, color = "#FFFFFFFF", linewidth = 0.5); Plot_Celltype_Sample
ggsave(plot = Plot_Celltype_Sample, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_Celltype_Sample.pdf", 
       width = 15, height = 10, dpi = 300, units = "in")

# Plot_Celltype_Group
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
Plot_Celltype_Group <- Seurat_Celltype_T_Cell@meta.data %>%
  dplyr::group_by(Group, Celltype_T_Cell) %>%
  dplyr::summarise(Count = n(), .groups = "drop_last") %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup() %>% 
  ggplot2::ggplot(mapping = aes(x = Group, y = Proportion, 
                                fill = Celltype_T_Cell, 
                                stratum = Celltype_T_Cell, 
                                alluvium = Celltype_T_Cell)) +
  ggplot2::scale_y_continuous(labels = scales::percent) +  
  ggplot2::scale_fill_manual(values = Colours_ScRNA_Cell) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Proportion of Cell Types (Group)",
                x = "Group", y = "Proportion", fill = "Cell Type") +
  ggplot2::theme(panel.grid = element_blank(), 
                 axis.title.x = element_text(size = 12, face = "bold"), 
                 axis.title.y = element_text(size = 12, face = "bold"), 
                 axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
                 axis.text.y = element_text(size = 10, face = "bold"), 
                 legend.text = element_text(size = 10, face = "bold"), 
                 legend.title = element_text(size = 12, face = "bold"), 
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                 legend.key.height = unit(3, "lines"),
                 legend.key.width  = unit(3, "lines")) + 
  ggplot2::geom_col(width = 0.5, color = NA, linewidth = 1) + 
  ggalluvial::geom_flow(width = 0.5, alpha = 0.25, knot.pos = 0.25, color = "#FFFFFFFF", linewidth = 0.5) + 
  ggalluvial::geom_alluvium(width = 0.5, alpha = 1, knot.pos = 0.25, fill = NA, color = "#FFFFFFFF", linewidth = 0.5); Plot_Celltype_Group
ggsave(plot = Plot_Celltype_Group, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_Celltype_Group.pdf", 
       width = 6, height = 10, dpi = 300, units = "in")

# Plot_Celltype
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
Plot_Celltype <- Seurat_Celltype_T_Cell@meta.data %>%
  dplyr::group_by(Celltype_T_Cell, Group) %>%
  dplyr::summarise(Count = n(), .groups = "drop_last") %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup() %>% 
  ggplot2::ggplot(mapping = aes(x = Celltype_T_Cell, y = Proportion, 
                                fill = Group, 
                                stratum = Group, 
                                alluvium = Group)) +
  ggplot2::scale_y_continuous(labels = scales::percent) +  
  ggplot2::scale_fill_manual(values = c("#FFA845FF", "#33CC33FF", "#1CCEFFFF")) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Proportion of Cell Types (Group)",
                x = "Group", y = "Proportion", fill = "Cell Type") +
  ggplot2::theme(panel.grid = element_blank(), 
                 axis.title.x = element_text(size = 12, face = "bold"), 
                 axis.title.y = element_text(size = 12, face = "bold"), 
                 axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
                 axis.text.y = element_text(size = 10, face = "bold"), 
                 legend.text = element_text(size = 10, face = "bold"), 
                 legend.title = element_text(size = 12, face = "bold"), 
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                 legend.key.height = unit(3, "lines"),
                 legend.key.width  = unit(3, "lines")) + 
  ggplot2::geom_col(width = 0.5, color = NA, linewidth = 1) + 
  ggalluvial::geom_flow(width = 0.5, alpha = 0.25, knot.pos = 0.25, color = "#FFFFFFFF", linewidth = 0.5) + 
  ggalluvial::geom_alluvium(width = 0.5, alpha = 1, knot.pos = 0.25, fill = NA, color = "#FFFFFFFF", linewidth = 0.5); Plot_Celltype
ggsave(plot = Plot_Celltype, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_Celltype.pdf", 
       width = 8, height = 10, dpi = 300, units = "in")

# Cell Ratio Radar ####
Plot_Radar_Group <- Seurat_Celltype_T_Cell@meta.data %>%
  dplyr::group_by(Group, Celltype_T_Cell) %>%
  dplyr::summarise(Count = n(), .groups = "drop_last") %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup() %>%
  dplyr::select(Group, Celltype_T_Cell, Proportion) %>%
  tidyr::pivot_wider(names_from = Celltype_T_Cell, values_from = Proportion, values_fill = 0) %>%
  dplyr::arrange(Group) %>% 
  ggradar2::ggradar2(base.size = 20,
                     # values.radar = c("0", "0.2", "0.4", "0.6", "0.8", "1"),
                     grid.line.width = 1, 
                     group.line.width = 1,
                     group.point.size = 2,
                     group.colours = c("#FFA845FF", "#33CC33FF", "#1CCEFFFF"),
                     fill = TRUE,
                     fill.alpha = 0.1,
                     axis.label.offset = 1.15,
                     axis.label.size = 4, 
                     # plot.title = "Proportion of Cell Types",
                     # legend.title = "Group",
                     legend.text.size = 12,
                     legend.position = "top",
                     background.circle.colour = "#FFFFFFFF"); Plot_Radar_Group
ggsave(plot = Plot_Radar_Group, 
       filename = "5_ScRNA_Celltype_T_Cell/Plot_Radar_Group.pdf", 
       width = 8, height = 8, dpi = 300, units = "in")

# Cell Ratio Boxplot
Plot_Cell_Box <- Seurat_Celltype_T_Cell@meta.data %>%
  dplyr::group_by(Sample, Celltype_T_Cell) %>%
  dplyr::summarise(Count = n(), .groups = "drop_last") %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(Group = dplyr::case_when(
    Sample %in% c("Nonresponse_1", "Nonresponse_2", "Nonresponse_3", "Nonresponse_4") ~ "Nonresponse",
    Sample %in% c("Response_1", "Response_2", "Response_3") ~ "Response",
    TRUE ~ "Unknown")) %>% 
  ggplot2::ggplot(mapping = aes(x = Group, y = Proportion, fill = Group)) + 
  ggplot2::stat_boxplot(geom = "errorbar", width = 0.5, linewidth = 1) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 1) +
  ggplot2::geom_jitter(mapping = aes(fill = Group), color = "#000000FF", shape = 21, size = 2, width = 0.2, alpha = 0.8) + 
  ggpubr::stat_compare_means(comparisons = list(c("Nonresponse", "Response")), 
                             method = "wilcox.test", 
                             label = "p.format", 
                             size = 5,
                             bracket.size = 1) +
  ggplot2::facet_wrap(~ Celltype_T_Cell, nrow = 1) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(title = "Celltype", y = "Cell ratio", x = "Group") +
  ggplot2::theme_classic(base_size = 10) + 
  ggplot2::theme(panel.grid = element_blank(), 
                 strip.background = element_rect(fill = "#FFFFFFFF", color = "#FFFFFFFF"), 
                 strip.text = element_text(size = 12, face = "bold"),
                 axis.title.x = element_text(size = 10, face = "bold"), 
                 axis.title.y = element_text(size = 10, face = "bold"), 
                 axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
                 axis.text.y = element_text(size = 10, face = "bold"), 
                 legend.text = element_text(size = 10, face = "bold"),
                 legend.title = element_text(size = 10, face = "bold"), 
                 plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
                 plot.subtitle = element_text(size = 15, face = "bold")) +
  ggplot2::scale_fill_manual(values = c("#FFA845FF", "#33CC33FF", "#1CCEFFFF")); Plot_Cell_Box
ggsave(plot = Plot_Cell_Box, filename = "5_ScRNA_Celltype_T_Cell/Plot_Cell_Box.pdf", 
       width = 10, height = 6, units = "in", dpi = 300)

# Cell Ratio Barplot
Plot_Cell_Bar <- Seurat_Celltype_T_Cell@meta.data %>%
  dplyr::group_by(Sample, Celltype_T_Cell) %>%
  dplyr::summarise(Count = n(), .groups = "drop_last") %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(Group = dplyr::case_when(
    Sample %in% c("Nonresponse_1", "Nonresponse_2", "Nonresponse_3", "Nonresponse_4") ~ "Nonresponse",
    Sample %in% c("Response_1", "Response_2", "Response_3") ~ "Response",
    TRUE ~ "Unknown")) %>% 
  ggplot2::ggplot(mapping = aes(x = Group,  y = Proportion, fill = Group)) + 
  ggplot2::stat_summary(fun = mean, geom = "bar",width = 0.6, color = "black", linewidth = 1, position = position_dodge(width = 0.6)) +
  ggplot2::stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, linewidth = 0.8, position = position_dodge(width = 0.6)) +
  ggplot2::geom_jitter(mapping = aes(fill = Group), color = "#000000FF", shape = 21, size = 2, width = 0.2, alpha = 0.8) + 
  ggplot2::facet_wrap(~ Celltype_T_Cell, nrow = 1) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(title = "Celltype", y = "Cell ratio", x = "Group") +
  ggplot2::theme_classic(base_size = 10) + 
  ggplot2::theme(panel.grid = element_blank(), 
                 strip.background = element_rect(fill = "#FFFFFFFF", color = "#FFFFFFFF"), 
                 strip.text = element_text(size = 12, face = "bold"),
                 axis.title.x = element_text(size = 10, face = "bold"), 
                 axis.title.y = element_text(size = 10, face = "bold"), 
                 axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
                 axis.text.y = element_text(size = 10, face = "bold"), 
                 legend.text = element_text(size = 10, face = "bold"),
                 legend.title = element_text(size = 10, face = "bold"), 
                 plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
                 plot.subtitle = element_text(size = 15, face = "bold")) +
  ggplot2::scale_fill_manual(values = c("#FFA845FF", "#33CC33FF", "#1CCEFFFF")); Plot_Cell_Bar
ggsave(plot = Plot_Cell_Bar, filename = "5_ScRNA_Celltype_T_Cell/Plot_Cell_Bar.pdf", 
       width = 10, height = 6, units = "in", dpi = 300)

# Log2FC [Nonresponse vs Response] ####
Log2FC_Nonresponse_Response <- Seurat_Celltype_T_Cell@meta.data %>%
  dplyr::group_by(Group, Celltype_T_Cell) %>%
  dplyr::summarise(Count = n(), .groups = "drop_last") %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Group %in% c("Nonresponse", "Response")) %>% 
  dplyr::select(Group, Celltype_T_Cell, Proportion) %>%
  tidyr::pivot_wider(names_from = Group, values_from = Proportion) %>%
  dplyr::mutate(Fold_Change = Nonresponse / Response,
                Log2FC = log2(Nonresponse / Response)) %>%
  ggplot2::ggplot(aes(x = Celltype_T_Cell, y = Log2FC, fill = Log2FC > 0)) +
  ggplot2::geom_bar(stat = "identity", width = 0.5) +
  ggplot2::geom_text(aes(label = round(Log2FC, 2),
                         vjust = ifelse(Log2FC > 0, -1, 2)), size = 3) +
  ggplot2::scale_fill_manual(values = c("TRUE" = "#FFA845FF", "FALSE" = "#33CC33FF"),
                             labels = c("Response > Nonresponse", "Nonresponse > Response")) +
  ggplot2::theme_bw(base_size = 15) +
  ggplot2::coord_cartesian(ylim = c(-3, 3)) +
  ggplot2::labs(title = "Nonresponse vs Response", x = "Celltype", y = "log2 Fold Change", fill = "log2 Fold Change") +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  ggplot2::theme(panel.grid = element_blank(), 
                 axis.title.x = element_text(size = 12, face = "bold"), 
                 axis.title.y = element_text(size = 12, face = "bold"), 
                 axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
                 axis.text.y = element_text(size = 10, face = "bold"), 
                 legend.text = element_text(size = 10, face = "bold"), 
                 legend.title = element_text(size = 12, face = "bold"), 
                 legend.position = "top",
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                 legend.key.size = unit(0.5, "in")); Log2FC_Nonresponse_Response
ggsave(plot = Log2FC_Nonresponse_Response, filename = "5_ScRNA_Celltype_T_Cell/Log2FC_Nonresponse_Response.pdf", 
       units = "in", width = 10, height = 6, dpi = 300)

# Add MetaData ####
load(file = "1_ScRNA_SeuratData/Seurat_Celltype_Major.RData")
head(Seurat_Celltype_Major@meta.data, n = 10); gc()
Seurat_Celltype_Major@meta.data <- Seurat_Celltype_Major@meta.data %>%
  tibble::rownames_to_column(var = "Cell_ID") %>%
  dplyr::left_join(Seurat_Celltype_T_Cell@meta.data %>%
                     dplyr::select(Celltype_T_Cell) %>% 
                     tibble::rownames_to_column(var = "Cell_ID"), by = "Cell_ID") %>%
  dplyr::mutate(Celltype_T_Cell = dplyr::case_when(
    !is.na(Celltype_T_Cell) ~ as.character(Celltype_T_Cell), 
    TRUE ~ as.character(Celltype_Major))) %>%
  tibble::column_to_rownames(var = "Cell_ID") %>%
  base::as.data.frame()
table(Seurat_Celltype_Major@meta.data$Celltype_Major)
table(Seurat_Celltype_Major@meta.data$Celltype_T_Cell)
save(Seurat_Celltype_Major, file = "1_ScRNA_SeuratData/Seurat_Celltype_Major.RData")
