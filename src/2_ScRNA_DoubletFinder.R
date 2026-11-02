# Set ####
options(timeout = 7355608); getOption("timeout")
options(stringsAsFactors = FALSE)
getwd(); rm(list = ls()); gc()
set.seed(19981108)
if (!dir.exists("2_ScRNA_DoubletFinder")) {dir.create("2_ScRNA_DoubletFinder")}

# Library ####
library(Seurat); packageVersion("Seurat")
library(tidyverse); packageVersion("tidyverse")
library(DoubletFinder); packageVersion("DoubletFinder")
library(crayon); packageVersion("crayon")
library(patchwork); packageVersion("patchwork")

# Load Data ####
load(file = "1_ScRNA_SeuratData/Seurat_List.RData")
Seurat_List; names(Seurat_List)

# DoubletFinder ####
Seurat_DoubletFinder_List <- list()
Time <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
for (Sample in names(Seurat_List)) {
  cat(red$bold(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Processing Samples: "), Sample, "|", Time(), "\n")
  Seurat <- Seurat_List[[Sample]]
  
  # NormalizeData & FindVariableFeatures & ScaleData & RunPCA & RunUMAP_TSNE & FindNeighbors & FindClusters
  Seurat <- Seurat %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(resolution = 1) %>% 
    RunUMAP(dims = 1:30) %>%
    RunTSNE(dims = 1:30)
  
  # Parameter Sweep
  Sweep_Res_List <- paramSweep(Seurat, PCs = 1:30, sct = FALSE)
  Sweep_Stats <- summarizeSweep(Sweep_Res_List, GT = FALSE)
  Bcmvn <- find.pK(Sweep_Stats)
  if (is.null(Bcmvn) || nrow(Bcmvn) == 0) {Best_Pk <- 0.09} else {Best_Pk <- Bcmvn$pK[which.max(Bcmvn$BCmetric)] %>%
      as.character() %>% as.numeric()}
  Homo_Typic_Prop <- modelHomotypic(Seurat$seurat_clusters)
  N_Exp_Poi <- round(0.075 * ncol(Seurat)) # Assuming 7.5% doublet formation rate - tailor for your dataset [0.075]
  N_Exp_Poi_Adj <- round(N_Exp_Poi * (1 - Homo_Typic_Prop))
  
  # DoubletFinder [reuse.pANN = FALSE]
  Seurat <- doubletFinder(seu = Seurat, 
                          PCs = 1:30, 
                          pN = 0.25, 
                          pK = Best_Pk, 
                          nExp = N_Exp_Poi_Adj, 
                          sct = FALSE)
  
  # Rename
  DF_Col <- grep("^DF.classifications", colnames(Seurat@meta.data), value = TRUE)
  pANN_Col <- grep("^pANN", colnames(Seurat@meta.data), value = TRUE)
  Seurat@meta.data <- Seurat@meta.data %>%
    dplyr::rename(DoubletFinder = dplyr::all_of(DF_Col),
                  pANN = dplyr::all_of(pANN_Col))
  # Return List
  Seurat_DoubletFinder_List[[Sample]] <- Seurat
  
  # Clear
  rm(Seurat, Sweep_Res_List, Sweep_Stats, Bcmvn, Best_Pk, Homo_Typic_Prop, N_Exp_Poi, N_Exp_Poi_Adj, DF_Col, pANN_Col)
  cat(red$bold(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Processing End: "), Sample, "|", Time(), "\n")
  gc()}
Seurat_DoubletFinder_List
save(Seurat_DoubletFinder_List, file = "2_ScRNA_DoubletFinder/Seurat_DoubletFinder_List.RData")

# UMAP DoubletFinder ####
Plot_UMAP_DoubletFinder_List <- list()
for (Sample in names(Seurat_DoubletFinder_List)) {
  Colours_ScRNA_DoubletFinder <- c("#FF6354FF", "#666666FF")
  Plot_UMAP_DoubletFinder <- DimPlot(object = Seurat_DoubletFinder_List[[Sample]],
                                     group.by = "DoubletFinder",
                                     pt.size = 0.1,
                                     raster = FALSE,
                                     reduction = "umap") +
    labs(title = Sample, 
         x = "UMAP 1", 
         y = "UMAP 2") + 
    theme_bw() +
    theme(legend.position = "top", 
          panel.grid = element_blank(),
          axis.title.x = element_text(size = 12, face = "bold"), 
          axis.title.y = element_text(size = 12, face = "bold"), 
          axis.text.x = element_text(size = 10, face = "bold"), 
          axis.text.y = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 10, face = "bold"), 
          legend.title = element_text(size = 12, face = "bold"), 
          plot.subtitle = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) + 
    scale_color_manual(values = Colours_ScRNA_DoubletFinder)
  Plot_UMAP_DoubletFinder_List[[Sample]] <- Plot_UMAP_DoubletFinder
  rm(Plot_UMAP_DoubletFinder, Sample, Colours_ScRNA_DoubletFinder)}
Plot_UMAP_DoubletFinder <- patchwork::wrap_plots(Plot_UMAP_DoubletFinder_List, ncol = 4)
Plot_UMAP_DoubletFinder
ggsave(plot = Plot_UMAP_DoubletFinder, 
       filename = "2_ScRNA_DoubletFinder/Plot_UMAP_DoubletFinder.pdf", 
       width = 16, height = 15, units = "in", dpi = 300) # [4:5]

# Merge ####
Seurat_DoubletFinder <- merge(x = Seurat_DoubletFinder_List[[1]],  
                              y = Seurat_DoubletFinder_List[-1], 
                              add.cell.ids = names(Seurat_DoubletFinder_List), 
                              project = "ScRNA") # [!]
head(Seurat_DoubletFinder@meta.data, n = 10); dim(Seurat_DoubletFinder); table(Seurat_DoubletFinder$DoubletFinder) # [Doublet: 8414 Singlet: 113959]
MetaData_DoubletFinder <- Seurat_DoubletFinder@meta.data %>% 
  dplyr::select(pANN, DoubletFinder)
head(MetaData_DoubletFinder, n = 10)
save(MetaData_DoubletFinder, file = "2_ScRNA_DoubletFinder/MetaData_DoubletFinder.RData")

# File Tree ####
fs::dir_tree("2_ScRNA_DoubletFinder/", max_depth = 2)

# End ####
sessionInfo()

