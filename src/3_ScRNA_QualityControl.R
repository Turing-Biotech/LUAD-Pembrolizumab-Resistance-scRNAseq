# Set ####
options(timeout = 7355608); getOption("timeout")
options(stringsAsFactors = FALSE)
getwd(); rm(list = ls()); gc()
set.seed(19981108)
if (!dir.exists("3_ScRNA_QualityControl")) {dir.create("3_ScRNA_QualityControl")}

# Library ####
library(Seurat); packageVersion("Seurat")
library(tidyverse); packageVersion("tidyverse")
library(scales); packageVersion("scales")
library(patchwork); packageVersion("patchwork")

# Load Data ####
load(file = "1_ScRNA_SeuratData/Seurat.RData")
head(Seurat@meta.data, n = 10)

# Genes Ratio ####
# Mitochondrial Genes
Percent_MT <- rownames(Seurat)[grep("^MT-", rownames(Seurat), ignore.case = TRUE)]
Percent_MT # [13]
Seurat <- PercentageFeatureSet(object = Seurat, 
                               features = Percent_MT, 
                               col.name = "Percent_MT")
fivenum(Seurat@meta.data$Percent_MT)
head(Seurat@meta.data, n = 10)

# Ribosomal Genes
Percent_RB <- rownames(Seurat)[grep("^RPL[0-9]+$|^RPS[0-9]+$", rownames(Seurat), ignore.case = TRUE)]
Percent_RB # [60]
Seurat <- PercentageFeatureSet(object = Seurat, 
                               features = Percent_RB,
                               col.name = "Percent_RB")
fivenum(Seurat@meta.data$Percent_RB)
head(Seurat@meta.data, n = 10)

# Hemoglobin Genes
Percent_HB <- rownames(Seurat)[grep("^HB[ABDGQZM][12]?$", rownames(Seurat), ignore.case = TRUE)]
Percent_HB # [9]
Seurat <- PercentageFeatureSet(object = Seurat, 
                               features = Percent_HB, 
                               col.name = "Percent_HB") 
fivenum(Seurat@meta.data$Percent_HB)
head(Seurat@meta.data, n = 10)

# [Before] Plot QC ####
table(Seurat@meta.data$Sample); gc()
Colours_ScRNA_Sample <- c(
  "#FF1C1CFF", "#FB8B05FF", "#FFFF00FF", "#33CC33FF", "#6699FFFF", "#1CFFD2FF", "#A246FFFF", "#FF7AC2FF")

# [Before] nFeature : nCount ####
Plot_nFeature_nCount <- FeatureScatter(object = Seurat, 
                                       feature1 = "nCount_RNA", 
                                       feature2 = "nFeature_RNA", 
                                       pt.size = 0.5, 
                                       group.by = "Sample", 
                                       raster = TRUE) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) + 
  scale_color_manual(values = Colours_ScRNA_Sample) + 
  scale_x_continuous(labels = scales::comma,
                     breaks = seq(0, 2000000, 10000)) + 
  scale_y_continuous(labels = scales::comma,
                     breaks = seq(0, 15000, 1000)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) +
  labs(x = "nCount RNA",
       y = "nFeature RNA")
Plot_nFeature_nCount

# [Before] Percent_MT : nCount ####
Plot_Percent_MT_nCount <- FeatureScatter(object = Seurat, 
                                         feature1 = "nCount_RNA", 
                                         feature2 = "Percent_MT", 
                                         pt.size = 0.5, 
                                         group.by = "Sample", 
                                         raster = TRUE) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) + 
  scale_color_manual(values = Colours_ScRNA_Sample) + 
  scale_x_continuous(labels = scales::comma, 
                     breaks = seq(0, 2000000, 10000)) + 
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 100, 10)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) +
  labs(x = "nCount RNA",
       y = "Percent MT")
Plot_Percent_MT_nCount

# [Before] nFeature ####
Plot_nFeature <- VlnPlot(object = Seurat, 
                         group.by = "Sample", 
                         features = "nFeature_RNA", 
                         assay = "RNA",
                         layer = "counts",
                         pt.size = 0, 
                         raster = TRUE) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        legend.key.size = unit(0.5, "in")) + 
  scale_fill_manual(values = Colours_ScRNA_Sample) + 
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 15000, 1000)) +
  labs(title = "nFeature RNA (Before Quality Control)",
       x = "Sample",
       y = "nFeature RNA")
Plot_nFeature

# [Before] nCount ####
Plot_nCount <- VlnPlot(object = Seurat, 
                       group.by = "Sample", 
                       features = "nCount_RNA", 
                       assay = "RNA",
                       layer = "counts",
                       pt.size = 0, 
                       raster = TRUE) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        legend.key.size = unit(0.5, "in")) + 
  scale_fill_manual(values = Colours_ScRNA_Sample) + 
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 2000000, 10000)) + 
  labs(title = "nCount RNA (Before Quality Control)",
       x = "Patient",
       y = "nCount RNA")
Plot_nCount

# [Before] Percent_MT ####
Plot_Percent_MT <- VlnPlot(object = Seurat, 
                           group.by = "Sample", 
                           features = "Percent_MT", 
                           assay = "RNA",
                           layer = "counts",
                           pt.size = 0, 
                           raster = TRUE) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        legend.key.size = unit(0.5, "in")) + 
  scale_fill_manual(values = Colours_ScRNA_Sample) + 
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 100, 5)) +
  labs(title = "Percent MT (Before Quality Control)",
       x = "Patient",
       y = "Percent MT")
Plot_Percent_MT

# [Before] Percent_HB ####
Plot_Percent_HB <- VlnPlot(object = Seurat, 
                           group.by = "Sample", 
                           features = "Percent_HB", 
                           assay = "RNA",
                           layer = "counts",
                           pt.size = 0, 
                           raster = TRUE) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        legend.key.size = unit(0.5, "in")) + 
  scale_fill_manual(values = Colours_ScRNA_Sample) + 
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 100, 5)) + 
  labs(title = "Percent HB (Before Quality Control)",
       x = "Patient",
       y = "Percent HB")
Plot_Percent_HB

# [Before] Percent_RB ####
Plot_Percent_RB <- VlnPlot(object = Seurat, 
                           group.by = "Sample", 
                           features = "Percent_RB", 
                           assay = "RNA",
                           layer = "counts",
                           pt.size = 0, 
                           raster = TRUE) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        legend.key.size = unit(0.5, "in")) + 
  scale_fill_manual(values = Colours_ScRNA_Sample) + 
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 100, 5)) + 
  labs(title = "Percent RB (Before Quality Control)",
       x = "Patient",
       y = "Percent RB")
Plot_Percent_RB

# [QC-Ing] Quality Control ####
load(file = "2_ScRNA_DoubletFinder/MetaData_DoubletFinder.RData")
dim(MetaData_DoubletFinder); dim(Seurat@meta.data); table(MetaData_DoubletFinder$DoubletFinder)

# DoubletFinder QC
MetaData_DoubletFinder <- MetaData_DoubletFinder[rownames(Seurat@meta.data), ]
Seurat@meta.data <- Seurat@meta.data %>% 
  tibble::rownames_to_column(var = "Barcode") %>% 
  dplyr::left_join(MetaData_DoubletFinder %>% 
                     tibble::rownames_to_column(var = "Barcode"), by = "Barcode") %>% 
  tibble::column_to_rownames(var = "Barcode") %>% 
  base::as.data.frame()
head(MetaData_DoubletFinder, n = 10); head(Seurat@meta.data, n = 10)

# Normal QC
Cells_QualityControl <- Seurat@meta.data %>% 
  dplyr::mutate(nFeature_nCount = dplyr::case_when(
    Sample == "Response_1" & nFeature_RNA < 4500 & nCount_RNA < 16000 ~ TRUE,
    Sample == "Response_2" & nFeature_RNA < 4500 & nCount_RNA < 15000 ~ TRUE,
    Sample == "Response_3" & nFeature_RNA < 2500 & nCount_RNA < 7000 ~ TRUE,
    Sample == "Nonresponse_1" & nFeature_RNA < 4500 & nCount_RNA < 15000 ~ TRUE,
    Sample == "Nonresponse_2" & nFeature_RNA < 5000 & nCount_RNA < 15000 ~ TRUE,
    Sample == "Nonresponse_3" & nFeature_RNA < 5000 & nCount_RNA < 20000 ~ TRUE,
    Sample == "Nonresponse_4" & nFeature_RNA < 5000 & nCount_RNA < 15000 ~ TRUE,
    TRUE ~ FALSE)) %>% 
  dplyr::filter(nFeature_nCount,
                DoubletFinder == "Singlet",
                Percent_MT < 10, 
                Percent_RB < 50,
                Percent_HB < 1,
                nFeature_RNA > 300,
                nCount_RNA > 500) %>%
  rownames()
Seurat_QualityControl <- subset(x = Seurat, cells = Cells_QualityControl)
print(dim(Seurat)) # 106385
print(dim(Seurat_QualityControl)) # 91576

# [After] nFeature : nCount ####
Plot_nFeature_nCount_QC <- FeatureScatter(object = Seurat_QualityControl, 
                                          feature1 = "nCount_RNA", 
                                          feature2 = "nFeature_RNA", 
                                          pt.size = 0.5, 
                                          group.by = "Sample", 
                                          raster = TRUE) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) + 
  scale_color_manual(values = Colours_ScRNA_Sample) + 
  scale_x_continuous(labels = scales::comma,
                     breaks = seq(0, 300000, 1000)) + 
  scale_y_continuous(labels = scales::comma,
                     breaks = seq(0, 15000, 500)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) +
  labs(x = "nCount RNA",
       y = "nFeature RNA")
Plot_nFeature_nCount_QC

# [After] Percent_MT : nCount ####
Plot_Percent_MT_nCount_QC <- FeatureScatter(object = Seurat_QualityControl, 
                                            feature1 = "nCount_RNA", 
                                            feature2 = "Percent_MT", 
                                            pt.size = 0.5, 
                                            group.by = "Sample", 
                                            raster = TRUE) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) + 
  scale_color_manual(values = Colours_ScRNA_Sample) + 
  scale_x_continuous(labels = scales::comma, 
                     breaks = seq(0, 200000, 1000)) + 
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 100, 1)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) +
  labs(x = "nCount RNA",
       y = "Percent MT")
Plot_Percent_MT_nCount_QC

# [After] nFeature ####
Plot_nFeature_QC <- VlnPlot(object = Seurat_QualityControl, 
                            group.by = "Sample", 
                            features = "nFeature_RNA", 
                            assay = "RNA",
                            layer = "counts",
                            pt.size = 0, 
                            raster = TRUE) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        legend.key.size = unit(0.5, "in")) + 
  scale_fill_manual(values = Colours_ScRNA_Sample) + 
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 15000, 500)) + 
  labs(title = "nFeature RNA (After Quality Control)",
       x = "Sample",
       y = "nFeature RNA")
Plot_nFeature_QC

# [After] nCount ####
Plot_nCount_QC <- VlnPlot(object = Seurat_QualityControl, 
                          group.by = "Sample", 
                          features = "nCount_RNA", 
                          assay = "RNA",
                          layer = "counts",
                          pt.size = 0, 
                          raster = TRUE) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        legend.key.size = unit(0.5, "in")) + 
  scale_fill_manual(values = Colours_ScRNA_Sample) + 
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 300000, 1000)) + 
  labs(title = "nCount RNA (After Quality Control)",
       x = "Patient",
       y = "nCount RNA")
Plot_nCount_QC

# [After] Percent_MT ####
Plot_Percent_MT_QC <- VlnPlot(object = Seurat_QualityControl, 
                              group.by = "Sample", 
                              features = "Percent_MT", 
                              assay = "RNA",
                              layer = "counts",
                              pt.size = 0, 
                              raster = TRUE) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        legend.key.size = unit(0.5, "in")) + 
  scale_fill_manual(values = Colours_ScRNA_Sample) + 
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 30, 1)) +
  labs(title = "Percent MT (After Quality Control)",
       x = "Patient",
       y = "Percent MT")
Plot_Percent_MT_QC

# [After] Percent_HB ####
Plot_Percent_HB_QC <- VlnPlot(object = Seurat_QualityControl, 
                              group.by = "Sample", 
                              features = "Percent_HB", 
                              assay = "RNA",
                              layer = "counts",
                              pt.size = 0, 
                              raster = TRUE) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        legend.key.size = unit(0.5, "in")) + 
  scale_fill_manual(values = Colours_ScRNA_Sample) + 
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 100, 5)) + 
  labs(title = "Percent HB (After Quality Control)",
       x = "Patient",
       y = "Percent HB")
Plot_Percent_HB_QC

# [After] Percent_RB ####
Plot_Percent_RB_QC <- VlnPlot(object = Seurat_QualityControl, 
                              group.by = "Sample", 
                              features = "Percent_RB", 
                              assay = "RNA",
                              layer = "counts",
                              pt.size = 0, 
                              raster = TRUE) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        legend.key.size = unit(0.5, "in")) + 
  scale_fill_manual(values = Colours_ScRNA_Sample) + 
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 50, 5)) + 
  labs(title = "Percent RB (After Quality Control)",
       x = "Patient",
       y = "Percent RB")
Plot_Percent_RB_QC

# Patchwork ####
Plot_QC_Scatter <- (Plot_nFeature_nCount + Plot_Percent_MT_nCount + 
                      Plot_nFeature_nCount_QC + Plot_Percent_MT_nCount_QC) + 
  plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(title = "Quality Control [nFeature vs nCount and Percent_MT vs nCount]",
                  theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
Plot_QC_Scatter
Plot_QC_Violin <- (Plot_nFeature + Plot_nCount + Plot_Percent_MT + Plot_Percent_HB + Plot_Percent_RB +
                     Plot_nFeature_QC + Plot_nCount_QC + Plot_Percent_MT_QC + Plot_Percent_HB_QC + Plot_Percent_RB_QC) + 
  plot_layout(ncol = 5, nrow = 2) + 
  plot_annotation(title = "Quality Control [nFeature, nCount, Percent_MT, Percent_HB, Percent_RB]",
                  theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
Plot_QC_Violin

# Save ####
ggsave(plot = Plot_QC_Scatter, 
       filename = "3_ScRNA_QualityControl/Plot_QC_Scatter.pdf",
       units = "in",
       width = 12,
       height = 10,
       dpi = 300)
ggsave(plot = Plot_QC_Violin, 
       filename = "3_ScRNA_QualityControl/Plot_QC_Violin.pdf",
       units = "in",
       limitsize = FALSE,
       width = 45,
       height = 15,
       dpi = 300)
save(Seurat_QualityControl, file = "1_ScRNA_SeuratData/Seurat_QualityControl.RData")
# load(file = "1_ScRNA_SeuratData/Seurat_QualityControl.RData")

# File Tree ####
fs::dir_tree("3_ScRNA_QualityControl/", max_depth = 2)

# End ####
sessionInfo()
