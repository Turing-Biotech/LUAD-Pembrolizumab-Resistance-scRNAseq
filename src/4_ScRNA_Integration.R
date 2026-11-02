# Set ####
options(timeout = 7355608); getOption("timeout")
options(stringsAsFactors = FALSE)
getwd(); rm(list = ls()); gc()
set.seed(19981108)
if (!dir.exists("4_ScRNA_Integration")) {dir.create("4_ScRNA_Integration")}

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

# Colours ####
Colours_ScRNA_Sample <- c(
  "#FF1C1CFF", "#FB8B05FF", "#FFFF00FF", "#33CC33FF", "#6699FFFF", "#1CFFD2FF", "#A246FFFF", "#FF7AC2FF")
Colours_ScRNA_Clustree <- c(
  "#FF6347FF", "#FB8B05FF", "#FFFF00FF", "#FCD337FF", "#66FF66FF", "#33CC33FF", 
  "#FF1C1CFF", "#1CCEFFFF", "#6699FFFF", "#1CFFD2FF", "#A246FFFF", "#FF7AC2FF")

# Load Data ####
load(file = "1_ScRNA_SeuratData/Seurat_QualityControl.RData")
Seurat_QualityControl

# NormalizeData ####
Seurat_Integration <- JoinLayers(object = Seurat_QualityControl) # JoinLayers
# load(file = "4_ScRNA_Integration/MetaData_Save.RData") # MetaData_Save
# Seurat_Integration <- subset(x = Seurat_Integration, cells = MetaData_Save) # MetaData_Save
Seurat_Integration <- NormalizeData(object = Seurat_Integration,
                                    normalization.method = "LogNormalize", 
                                    scale.factor = 10000,
                                    verbose = TRUE)
Seurat_Integration; rm(Seurat_QualityControl); gc(); dim(Seurat_Integration)

# FindVariableFeatures ####
Seurat_Integration <- FindVariableFeatures(object = Seurat_Integration, 
                                           selection.method = "vst", 
                                           nfeatures = 3000, # [!]
                                           verbose = TRUE)
HVGs_Integration <- VariableFeatures(object = Seurat_Integration)
HVGs_Integration

# Plot_HVGs ####
Plot_HVGs <- VariableFeaturePlot(object = Seurat_Integration) + 
  NoLegend() + 
  theme_classic()+ 
  theme(panel.grid = element_blank(),
        legend.position = "none")
Plot_HVGs_Top <- LabelPoints(plot = Plot_HVGs, 
                             points = head(VariableFeatures(Seurat_Integration), 20), 
                             repel = TRUE) + 
  theme_classic()+ 
  theme(panel.grid = element_blank(), 
        axis.title.y = element_blank())
Plot_HVGs_Top
Plot_HVGs_All <- (Plot_HVGs + Plot_HVGs_Top) + 
  plot_layout(ncol = 2, nrow = 1) +  
  plot_annotation(title = "HVGs", theme = theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold")))
Plot_HVGs_All
ggsave(plot = Plot_HVGs_All, 
       filename = "4_ScRNA_Integration/Plot_HVGs_All.pdf", 
       width = 8, height = 8, dpi = 300, units = "in")

# CellCycle ####
S_Genes_Human <- Seurat::cc.genes.updated.2019$s.genes
G2M_Genes_Human <- Seurat::cc.genes.updated.2019$g2m.genes

# CaseMatch
S_Genes_Match <- CaseMatch(S_Genes_Human, rownames(Seurat_Integration))
G2M_Genes_Match <- CaseMatch(G2M_Genes_Human, rownames(Seurat_Integration))

Seurat_Integration <- CellCycleScoring(object = Seurat_Integration,
                                       s.features = S_Genes_Match,
                                       g2m.features = G2M_Genes_Match,
                                       set.ident = FALSE)
head(Seurat_Integration@meta.data, n = 10)
table(Seurat_Integration@meta.data$Phase)

# ScaleData ####
Seurat_Integration <- ScaleData(object = Seurat_Integration, 
                                vars.to.regress = c("S.Score", "G2M.Score", "Percent_MT"),
                                features = VariableFeatures(object = Seurat_Integration))
Seurat_Integration
DefaultAssay(Seurat_Integration)

# SCT ####
Seurat_Integration <- SCTransform(object = Seurat_Integration, 
                                  assay = "RNA",
                                  new.assay.name = "SCT",
                                  variable.features.n = 3000, # [!]             
                                  vars.to.regress = c("S.Score", "G2M.Score", "Percent_MT"),
                                  verbose = TRUE)
DefaultAssay(Seurat_Integration)

# PCA ####
Seurat_Integration <- RunPCA(object = Seurat_Integration, 
                             reduction.name = "pca",
                             features = VariableFeatures(object = Seurat_Integration),
                             npcs = 50, 
                             seed.use = 19981108,
                             verbose = TRUE)
names(Seurat_Integration@reductions)
ElbowPlot(Seurat_Integration, ndims = 50, reduction = "pca")
Seurat_Integration

# Dims
PC_Stdev <- Seurat_Integration@reductions$pca@stdev
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
names(Seurat_Integration@reductions)
Seurat_Integration <- FindNeighbors(object = Seurat_Integration, 
                                    reduction = "pca", 
                                    dims = Dims)
Seurat_Integration <- FindClusters(object = Seurat_Integration, 
                                   resolution = 1, 
                                   cluster.name = "PCA_1",
                                   random.seed = 19981108,
                                   algorithm = 4)
head(Seurat_Integration@meta.data, n = 10)

# UMAP & TSNE [PCA]
# UMAP
Seurat_Integration <- RunUMAP(object = Seurat_Integration, 
                              dims = Dims,
                              min.dist = 0.3, # 0.3 [-]
                              n.neighbors = 50L, # 30L [+]
                              spread = 0.6, # 1 [-]
                              reduction.name = "pca_umap",
                              reduction = "pca")
names(Seurat_Integration@reductions)

# Plot_PCA_Sample
Plot_PCA_Sample <- DimPlot(object = Seurat_Integration, 
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
  guides(color = guide_legend(override.aes = list(size = 8))) +
  scale_color_manual(values = Colours_ScRNA_Sample)
Plot_PCA_Sample

# Plot_UMAP:PCA
Plot_UMAP_PCA <- DimPlot(object = Seurat_Integration, 
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
  guides(color = guide_legend(override.aes = list(size = 8))) +
  scale_color_manual(values = Colours_ScRNA_Sample)
Plot_UMAP_PCA

# Harmony ####
DefaultAssay(Seurat_Integration)
names(Seurat_Integration@reductions)
Seurat_Integration <- RunHarmony(Seurat_Integration,
                                 group.by.vars = "orig.ident",
                                 reduction.use = "pca",
                                 dims = Dims, 
                                 reduction.save = "harmony",
                                 project.dim = FALSE,
                                 plot_convergence = TRUE)

# FindNeighbors & FindClusters
Seurat_Integration <- FindNeighbors(object = Seurat_Integration, 
                                    reduction = "harmony", 
                                    dims = Dims)
Seurat_Integration <- FindClusters(object = Seurat_Integration, 
                                   resolution = 1, 
                                   cluster.name = "Harmony_1",
                                   random.seed = 19981108,
                                   algorithm = 4)
print(names(Seurat_Integration@reductions))

# UMAP & TSNE [Harmony]
Seurat_Integration <- RunUMAP(object = Seurat_Integration, 
                              dims = Dims,
                              min.dist = 0.3, # 0.3 [-]
                              n.neighbors = 50L, # 30L [+]
                              spread = 0.6, # 1 [-]
                              reduction.name = "harmony_umap",
                              reduction = "harmony")
print(names(Seurat_Integration@reductions))

# Plot_Harmony_Sample
Plot_Harmony_Sample <- DimPlot(object = Seurat_Integration, 
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
  guides(color = guide_legend(override.aes = list(size = 8))) +
  scale_color_manual(values = Colours_ScRNA_Sample)
Plot_Harmony_Sample

# Plot_UMAP:Harmony
Plot_UMAP_Harmony <- DimPlot(object = Seurat_Integration, 
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
  guides(color = guide_legend(override.aes = list(size = 8))) +
  scale_color_manual(values = Colours_ScRNA_Sample)
Plot_UMAP_Harmony

# Patchwork ####
Plot_UMAP_PCA <- Plot_PCA_Sample | Plot_UMAP_PCA; Plot_UMAP_PCA
Plot_UMAP_Harmony <- Plot_Harmony_Sample | Plot_UMAP_Harmony; Plot_UMAP_Harmony
ggsave(plot = Plot_UMAP_PCA, 
       filename = "4_ScRNA_Integration/Plot_UMAP_PCA.pdf",
       width = 10, height = 4, units = "in", dpi = 300)
ggsave(plot = Plot_UMAP_Harmony, 
       filename = "4_ScRNA_Integration/Plot_UMAP_Harmony.pdf",
       width = 10, height = 4, units = "in", dpi = 300)

# FindNeighbors & FindClusters ####
names(Seurat_Integration@reductions)
# Seurat_Integration <- FindNeighbors(object = Seurat_Integration, reduction = "harmony", dims = Dims, k.param = 20) # Harmony
Seurat_Integration <- FindNeighbors(object = Seurat_Integration, reduction = "harmony", dims = Dims, k.param = 10) # Harmony

# FindClusters ####
Resolution <- c(0.01, 0.05, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)
for (Res in Resolution) {
  Seurat_Integration <- FindClusters(object = Seurat_Integration, 
                                     resolution = Res, 
                                     cluster.name = paste0("Resolution_", Res),
                                     random.seed = 19981108,
                                     algorithm = 4)}
head(Seurat_Integration@meta.data, n = 10)
for (Res in Resolution) {
  colname <- paste0("Resolution_", Res)
  cat("Resolution:", Res, "\n")
  print(table(Seurat_Integration[[colname]]))}

# UMAP ####
Seurat_Integration <- RunUMAP(object = Seurat_Integration, 
                              dims = Dims, 
                              min.dist = 0.3, # 0.3 [-]
                              n.neighbors = 50L, # 30L [+]
                              spread = 0.6, # 1 [-]
                              reduction.name = "harmony_umap_cluster", # !!!
                              reduction = "harmony") # !!!
names(Seurat_Integration@reductions)

# Plot UMAP Resolution ####
names(Seurat_Integration@reductions)
# Plot_UMAP_Resolution
Plot_UMAP_Resolution <- cowplot::plot_grid(
  DimPlot(object = Seurat_Integration, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_0.01") + ggtitle("Harmony: 0.01"),
  DimPlot(object = Seurat_Integration, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_0.05") + ggtitle("Harmony: 0.05"),
  DimPlot(object = Seurat_Integration, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_0.2") + ggtitle("Harmony: 0.2"),
  DimPlot(object = Seurat_Integration, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_0.4") + ggtitle("Harmony: 0.4"),
  DimPlot(object = Seurat_Integration, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_0.6") + ggtitle("Harmony: 0.6"),
  DimPlot(object = Seurat_Integration, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_0.8") + ggtitle("Harmony: 0.8"),
  DimPlot(object = Seurat_Integration, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_1")   + ggtitle("Harmony: 1.0"),
  DimPlot(object = Seurat_Integration, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_1.2") + ggtitle("Harmony: 1.2"),
  DimPlot(object = Seurat_Integration, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_1.4") + ggtitle("Harmony: 1.4"),
  DimPlot(object = Seurat_Integration, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_1.6") + ggtitle("Harmony: 1.6"),
  DimPlot(object = Seurat_Integration, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_1.8") + ggtitle("Harmony: 1.8"),
  DimPlot(object = Seurat_Integration, raster = TRUE, reduction = "harmony_umap_cluster", label = TRUE, group.by = "Resolution_2") + ggtitle("Harmony: 2.0"),
  ncol = 4)
ggsave(plot = Plot_UMAP_Resolution, 
       filename = "4_ScRNA_Integration/Plot_UMAP_Resolution.pdf",
       width = 20, height = 12, units = "in", dpi = 300)

# Plot_Clustree
head(Seurat_Integration@meta.data, n = 10)
Prefix <- Seurat_Integration@meta.data %>% 
  tibble::as_tibble(rownames = "Cell_ID") %>% 
  dplyr::select(matches("^Resolution_[0-9.]+$")) %>% 
  names()
Prefix
Plot_Clustree <- Seurat_Integration@meta.data %>%
  dplyr::select(all_of(Prefix)) %>%
  clustree::clustree(prefix = "Resolution_") + 
  labs(colour = "Harmony") + 
  guides(color = guide_legend(override.aes = list(size = 8))) +
  scale_color_manual(values = Colours_ScRNA_Clustree)
Plot_Clustree
ggsave(plot = Plot_Clustree, 
       filename = "4_ScRNA_Integration/Plot_Clustree.pdf", 
       width = 15, height = 10, units = "in", dpi = 300)

# Resolution #### 
Seurat_Integration@meta.data <- Seurat_Integration@meta.data %>%
  dplyr::mutate(dplyr::across(all_of(grep("^Resolution_", colnames(Seurat_Integration@meta.data), value = TRUE)),
                              ~ factor(as.integer(as.character(.)), levels = sort(unique(as.integer(as.character(.)))))))
str(Seurat_Integration@meta.data)
for (Res in Resolution) {
  colname <- paste0("Resolution_", Res)
  cat("Resolution:", Res, "\n")
  print(table(Seurat_Integration[[colname]]))}

# UMAP & TSNE Coord ####
head(Seurat_Integration@meta.data, n = 10)
class(Seurat_Integration@meta.data)
head(Seurat_Integration@reductions$harmony_umap_cluster@cell.embeddings, n = 10) # harmony_umap_cluster
class(Seurat_Integration@reductions$harmony_umap_cluster@cell.embeddings)

# Extract
UMAP_Coord <- Seurat_Integration@reductions$harmony_umap_cluster@cell.embeddings %>% 
  base::as.data.frame() %>% 
  dplyr::rename(UMAP_1 = harmonyumapcluster_1) %>% 
  dplyr::rename(UMAP_2 = harmonyumapcluster_2)
all(rownames(Seurat_Integration@meta.data)  == rownames(UMAP_Coord))
which(rownames(Seurat_Integration@meta.data) != rownames(UMAP_Coord))

# Bind
Seurat_Integration@meta.data <- Seurat_Integration@meta.data %>% 
  dplyr::bind_cols(UMAP_Coord) 
head(Seurat_Integration@meta.data, n = 10)
class(Seurat_Integration@meta.data)

# Save ####
head(Seurat_Integration@meta.data, n = 10)
save(Seurat_Integration, file = "1_ScRNA_SeuratData/Seurat_Integration.RData")
# load(file = "1_ScRNA_SeuratData/Seurat_Integration.RData")

# File Tree ####
fs::dir_tree("4_ScRNA_Integration/", max_depth = 2)

# End ####
sessionInfo()


