# Set ####
options(timeout = 7355608); getOption("timeout")
options(stringsAsFactors = FALSE)
getwd(); rm(list = ls()); gc()
set.seed(19981108)
if (!dir.exists("5_ScRNA_Celltype_Major")) {dir.create("5_ScRNA_Celltype_Major")}

# Library ####
library(Seurat); packageVersion("Seurat")
library(tidyverse); packageVersion("tidyverse")
library(patchwork); packageVersion("patchwork")
library(SingleR); packageVersion("SingleR")
library(celldex); packageVersion("celldex")
library(ggrepel); packageVersion("ggrepel")
library(ggrastr); packageVersion("ggrastr")
library(ggalluvial); packageVersion("ggalluvial")
library(ggradar2); packageVersion("ggradar2")
library(ggpubr); packageVersion("ggpubr")
library(cols4all); packageVersion("cols4all")

# Colours_Cell ####
Colours_ScRNA_Sample <- c(
  "#FF1C1CFF", "#FB8B05FF", "#FFFF00FF", "#33CC33FF", "#6699FFFF", "#1CFFD2FF", "#A246FFFF", "#FF7AC2FF")
Colours_ScRNA_Cell <- c(
  "#33CC33FF", "#1CFFD2FF", "#1CCEFFFF", "#FFA845FF", "#C4A1FFFF", 
  "#FF6354FF", "#FF7AC2FF", "#9933CCFF", "#FFEB4EFF", "#FFC6FFFF")
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
load(file = "1_ScRNA_SeuratData/Seurat_Integration.RData")
Seurat_Celltype_Major <- Seurat_Integration; rm(Seurat_Integration); gc()
# view(Seurat_Celltype_Major@meta.data)

# Resolution ####
Seurat_Celltype_Major
DefaultAssay(Seurat_Celltype_Major)
names(Seurat_Celltype_Major@reductions)
names(Seurat_Celltype_Major@assays$RNA@layers)
table(Seurat_Celltype_Major@active.ident) 
head(Seurat_Celltype_Major@meta.data, n = 10) 
Resolution_Select <- "Resolution_1" # !!!
Idents(Seurat_Celltype_Major) <- Resolution_Select
table(Seurat_Celltype_Major@active.ident) 
names(Seurat_Celltype_Major@reductions)

# UMAP & TSNE Cluster ####
source(file = "0_Axes.R")
MetaData_Cluster <- Seurat_Celltype_Major@meta.data
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
       filename = "5_ScRNA_Celltype_Major/Plot_UMAP_Cluster.pdf", 
       width = 9, height = 6, dpi = 300, units = "in")

# FindAllMarkers [Cluster] ####
head(Seurat_Celltype_Major@meta.data, n = 10)
FindAllMarkers_Cluster <- FindAllMarkers(object = Seurat_Celltype_Major, 
                                         only.pos = TRUE,
                                         group.by = Resolution_Select,
                                         min.pct = 0.2,
                                         test.use = "wilcox", 
                                         logfc.threshold = 0.5,
                                         verbose = TRUE)
write_csv(x = FindAllMarkers_Cluster, file = "5_ScRNA_Celltype_Major/FindAllMarkers_Cluster.csv")

# Marker Cheak ####
Markers_Cheak <- tibble::tribble(
  ~Celltype, ~Markers,
  "B_Cell", "CD79A", 
  "B_Cell", "CD79B", 
  "B_Cell", "BANK1", 
  "B_Cell", "MS4A1", 
  
  "T_Cell", "CD3D", 
  "T_Cell", "CD3E", 
  "T_Cell", "CD3G", 
  "T_Cell", "TRAC", 
  "T_Cell", "TRBC1", 
  "T_Cell", "TRBC2", 
  "T_Cell", "CD247", 
  "T_Cell", "CCR7", 
  "T_Cell", "FOXP3", 
  
  "NK_Cell", "NKG7", 
  "NK_Cell", "NCAM1", 
  "NK_Cell", "KLRD1", 
  "NK_Cell", "KLRF1",
  "NK_Cell", "GNLY", 
  "NK_Cell", "GZMB",
  
  "Platelet", "PPBP", 
  "Platelet", "PF4", 
  "Platelet", "ITGA2B", 
  "Platelet", "ITGB3", 
  "Platelet", "GP9", 
  "Platelet", "SPARC", 
  
  "Plasma_Cell", "JCHAIN", # [CD138] Joining Chain Of Multimeric IgA And IgM
  "Plasma_Cell", "IGHA1", # Immunoglobulin Heavy Constant Alpha 1
  "Plasma_Cell", "MZB1", # Marginal Zone B And B1 Cell Specific Protein
  "Plasma_Cell", "DERL3", # Derlin 3

  "Monocyte", "CD14", 
  "Monocyte", "LYZ", 
  "Monocyte", "FCN1", 
  "Monocyte", "S100A8", 
  "Monocyte", "S100A9", 
  
  # End 
  )

# Markers Cheak [Plot_Dot_Markers_Cheak] ####
Plot_Dot_Markers_Cheak <- ggplot(data = DotPlot(object = Seurat_Celltype_Major, 
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
       filename = "5_ScRNA_Celltype_Major/Plot_Dot_Markers_Cheak.pdf", 
       width = 25, height = 15, dpi = 300, units = "in")

# Cell Save ####
# MetaData_Save <- Seurat_Celltype_Major@meta.data %>% 
#   dplyr::filter(!(Resolution_1 %in% c("", ""))) %>% 
#   base::rownames()
# MetaData_Save
# save(MetaData_Save, file = "4_ScRNA_Integration/MetaData_Save.RData")

# Marker Annotation ####
# Resolution_1 %in% c("", "") ~ "", 
Seurat_Celltype_Major@meta.data <- Seurat_Celltype_Major@meta.data %>%
  dplyr::mutate(Celltype_Major = dplyr::case_when(
    Resolution_1 %in% c("8", "24") ~ "B_Cell", 
    Resolution_1 %in% c("28") ~ "Plasma_Cell", 
    Resolution_1 %in% c("1", "2", "3", "6", "12", "14", "16", "22", "29") ~ "T_Cell", 
    Resolution_1 %in% c("9", "20") ~ "NK_Cell", 
    Resolution_1 %in% c("10") ~ "Platelet", 
    Resolution_1 %in% c("4", "5", "7", "11", "13", "15",  "17", "18", "19", "21", "23", "25", "26", "27") ~ "Monocyte", 
    TRUE ~ as.character(Resolution_1)))
head(Seurat_Celltype_Major@meta.data, n = 10)
table(Seurat_Celltype_Major@meta.data$Celltype_Major)

# Test ####
Plot_UMAP_TEST <- DimPlot(object = Seurat_Celltype_Major, 
                          reduction = "harmony_umap_cluster", # !!!
                          group.by = "Celltype_Major",  
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
       filename = "5_ScRNA_Celltype_Major/Plot_UMAP_TEST.pdf", 
       width = 8, height = 6, dpi = 300, units = "in")

# Optimization ####
# UMAP_1 >= X & UMAP_1 <= X & UMAP_2 >= X & UMAP_2 <= X & Celltype_Major == "X" ~ "X", 
Seurat_Celltype_Major@meta.data <- Seurat_Celltype_Major@meta.data %>%
  dplyr::mutate(Celltype_Major = dplyr::case_when(
UMAP_1 >= X & UMAP_1 <= X & UMAP_2 >= X & UMAP_2 <= X & Celltype_Major == "X" ~ "X", 
    TRUE ~ Celltype_Major))

# Marker Celltype ####
Markers_Celltype <- tibble::tribble(
  ~Celltype, ~Markers,
  "B_Cell", "CD79A", 
  "B_Cell", "CD79B", 
  "B_Cell", "BANK1", 
  "B_Cell", "MS4A1", 
  
  "T_Cell", "CD3D", 
  "T_Cell", "CD3E", 
  "T_Cell", "CD3G", 
  
  "NK_Cell", "NCAM1", 
  "NK_Cell", "KLRD1", 
  "NK_Cell", "KLRF1",
  
  "Platelet", "PPBP", 
  "Platelet", "PF4", 
  "Platelet", "SPARC", 
  
  "Plasma_Cell", "JCHAIN", # [CD138] Joining Chain Of Multimeric IgA And IgM
  "Plasma_Cell", "IGHA1", # Immunoglobulin Heavy Constant Alpha 1
  "Plasma_Cell", "MZB1", # Marginal Zone B And B1 Cell Specific Protein

  "Monocyte", "CD14", 
  "Monocyte", "LYZ", 
  "Monocyte", "FCN1", 
  "Monocyte", "S100A8", 
  "Monocyte", "S100A9", 
  )

# Markers Celltype [Plot_Dot_Markers_Celltype] ####
Plot_Dot_Markers_Celltype <- ggplot(data = DotPlot(object = Seurat_Celltype_Major, 
                                                   group.by = "Celltype_Major", 
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
       filename = "5_ScRNA_Celltype_Major/Plot_Dot_Markers_Celltype.pdf", 
       width = 10, height = 5, dpi = 300, units = "in")

# UMAP & TSNE Marker ####
source("0_Axes.R")
MetaData_Marker <- Seurat_Celltype_Major@meta.data
MetaData_Marker$MS4A1 <- FetchData(Seurat_Celltype_Major, vars = "MS4A1")[, 1] # B_Cell
MetaData_Marker$JCHAIN <- FetchData(Seurat_Celltype_Major, vars = "JCHAIN")[, 1] # Plasma_Cell
MetaData_Marker$PPBP <- FetchData(Seurat_Celltype_Major, vars = "PPBP")[, 1] # Platelet
MetaData_Marker$CD3D <- FetchData(Seurat_Celltype_Major, vars = "CD3D")[, 1] # T_Cell
MetaData_Marker$KLRF1 <- FetchData(Seurat_Celltype_Major, vars = "KLRF1")[, 1] # NK_Cell
MetaData_Marker$LYZ <- FetchData(Seurat_Celltype_Major, vars = "LYZ")[, 1] # Monocyte
head(MetaData_Marker, n = 10)
MetaData_Marker <- MetaData_Marker %>% 
  dplyr::select(UMAP_1, UMAP_2,  
                MS4A1, JCHAIN, PPBP, CD3D, KLRF1, LYZ) %>% 
  tibble::as_tibble() %>% 
  tidyr::pivot_longer(cols = -c(UMAP_1, UMAP_2),
                      names_to = "Marker", values_to = "expression")
head(MetaData_Marker, n = 10)

# UMAP
Plot_UMAP_Marker <- ggplot(data = MetaData_Marker,
                           mapping = aes(x = UMAP_1, y = UMAP_2, color = expression)) +
  geom_point_rast(size = 0, alpha = 0.8, raster.dpi = 300, scale = 0.5) +
  scale_color_gradientn(colors = c("#D9D9D9FF", "#FF1C1CFF"), 
                        name = "Expression") +
  labs(title = "Marker Expression",x = "", y = "") +
  facet_wrap(~ Marker , ncol = 3) + 
  NoAxes() +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_rect(fill = "#FFFFFFFF", color = "#FFFFFFFF"), 
        strip.text = element_text(size = 12, face = "bold"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
Plot_UMAP_Marker <- Add_Left_Axes(Plot = Plot_UMAP_Marker, Meta_Data = MetaData_Marker,
                                  X_Col = "UMAP_1", Y_Col = "UMAP_2",
                                  X_Label = "UMAP 1", Y_Label = "UMAP 2",
                                  X_Label_Offset = 1, Y_Label_Offset = 1); Plot_UMAP_Marker
ggsave(plot = Plot_UMAP_Marker, 
       filename = "5_ScRNA_Celltype_Major/Plot_UMAP_Marker.pdf", 
       width = 13, height = 8, dpi = 300, units = "in")

# UMAP & TSNE Celltype ####
source(file = "0_Axes.R")
MetaData_Celltype_Major <- Seurat_Celltype_Major@meta.data
MetaData_Celltype_Major
save(MetaData_Celltype_Major, file = "5_ScRNA_Celltype_Major/MetaData_Celltype_Major.RData")

# UMAP Celltype
Plot_UMAP_Celltype <- ggplot(data = MetaData_Celltype_Major %>% # !!!
                               dplyr::group_by(Celltype_Major) %>% # !!!
                               dplyr::mutate(Label = paste0(Celltype_Major, " (", n(), ")")) %>% # !!!
                               dplyr::mutate(Label = as.factor(Label)) %>% 
                               dplyr::ungroup(), 
                             mapping = aes(x = UMAP_1, # !!!
                                           y = UMAP_2, # !!!
                                           color = Label)) +
  geom_point_rast(size = 0, alpha = 0.8, raster.dpi = 300, scale = 0.5) + 
  geom_text_repel(data = MetaData_Celltype_Major %>%
                    dplyr::group_by(Label = as.factor(Celltype_Major)) %>% # !!!
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
  labs(title = "UMAP (Celltype_Major)", x = "", y = "", color = "Celltype_Major") + # !!!
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
Plot_UMAP_Celltype <- Add_Left_Axes(Plot = Plot_UMAP_Celltype, Meta_Data = MetaData_Celltype_Major, 
                                    X_Col = "UMAP_1", Y_Col = "UMAP_2", 
                                    X_Label = "UMAP 1", Y_Label = "UMAP 2", 
                                    X_Label_Offset = 1, Y_Label_Offset = 1); Plot_UMAP_Celltype
ggsave(plot = Plot_UMAP_Celltype, 
       filename = "5_ScRNA_Celltype_Major/Plot_UMAP_Celltype.pdf", 
       width = 9, height = 6, dpi = 300, units = "in")

# FindAllMarkers [Celltype] ####
head(Seurat_Celltype_Major@meta.data, n = 10)
FindAllMarkers_Celltype <- FindAllMarkers(object = Seurat_Celltype_Major, 
                                          only.pos = TRUE,
                                          group.by = "Celltype_Major",
                                          min.pct = 0.2,
                                          test.use = "wilcox", 
                                          logfc.threshold = 0.5,
                                          verbose = TRUE)
write_csv(x = FindAllMarkers_Celltype, file = "5_ScRNA_Celltype_Major/FindAllMarkers_Celltype.csv")

# Cell Ratio Stacked Bar Plot ####
head(Seurat_Celltype_Major@meta.data, n = 10)
Plot_Celltype_Sample <- Seurat_Celltype_Major@meta.data %>%
  dplyr::group_by(Sample, Celltype_Major) %>%
  dplyr::summarise(Count = n(), .groups = "drop_last") %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup() %>% 
  ggplot2::ggplot(mapping = aes(x = Sample, y = Proportion, 
                                fill = Celltype_Major,
                                stratum = Celltype_Major, 
                                alluvium = Celltype_Major)) +
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
       filename = "5_ScRNA_Celltype_Major/Plot_Celltype_Sample.pdf", 
       width = 15, height = 10, dpi = 300, units = "in")

# Plot_Celltype_Group
head(Seurat_Celltype_Major@meta.data, n = 10)
Plot_Celltype_Group <- Seurat_Celltype_Major@meta.data %>%
  dplyr::group_by(Group, Celltype_Major) %>%
  dplyr::summarise(Count = n(), .groups = "drop_last") %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup() %>% 
  ggplot2::ggplot(mapping = aes(x = Group, y = Proportion, 
                                fill = Celltype_Major, 
                                stratum = Celltype_Major, 
                                alluvium = Celltype_Major)) +
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
       filename = "5_ScRNA_Celltype_Major/Plot_Celltype_Group.pdf", 
       width = 6, height = 10, dpi = 300, units = "in")

# Plot_Celltype
head(Seurat_Celltype_Major@meta.data, n = 10)
Plot_Celltype <- Seurat_Celltype_Major@meta.data %>%
  dplyr::group_by(Celltype_Major, Group) %>%
  dplyr::summarise(Count = n(), .groups = "drop_last") %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup() %>% 
  ggplot2::ggplot(mapping = aes(x = Celltype_Major, y = Proportion, 
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
       filename = "5_ScRNA_Celltype_Major/Plot_Celltype.pdf", 
       width = 8, height = 10, dpi = 300, units = "in")

# Cell Ratio Radar ####
Plot_Radar_Group <- Seurat_Celltype_Major@meta.data %>%
  dplyr::group_by(Group, Celltype_Major) %>%
  dplyr::summarise(Count = n(), .groups = "drop_last") %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup() %>%
  dplyr::select(Group, Celltype_Major, Proportion) %>%
  tidyr::pivot_wider(names_from = Celltype_Major, values_from = Proportion, values_fill = 0) %>%
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
       filename = "5_ScRNA_Celltype_Major/Plot_Radar_Group.pdf", 
       width = 8, height = 8, dpi = 300, units = "in")

# Cell Ratio Boxplot
Plot_Cell_Box <- Seurat_Celltype_Major@meta.data %>%
  dplyr::group_by(Sample, Celltype_Major) %>%
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
  ggplot2::facet_wrap(~ Celltype_Major, nrow = 1) +
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
ggsave(plot = Plot_Cell_Box, filename = "5_ScRNA_Celltype_Major/Plot_Cell_Box.pdf", 
       width = 10, height = 6, units = "in", dpi = 300)

# Cell Ratio Barplot
Plot_Cell_Bar <- Seurat_Celltype_Major@meta.data %>%
  dplyr::group_by(Sample, Celltype_Major) %>%
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
  ggplot2::facet_wrap(~ Celltype_Major, nrow = 1) +
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
ggsave(plot = Plot_Cell_Bar, filename = "5_ScRNA_Celltype_Major/Plot_Cell_Bar.pdf", 
       width = 10, height = 6, units = "in", dpi = 300)

# Log2FC [Nonresponse vs Response] ####
Log2FC_Nonresponse_Response <- Seurat_Celltype_Major@meta.data %>%
  dplyr::group_by(Group, Celltype_Major) %>%
  dplyr::summarise(Count = n(), .groups = "drop_last") %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Group %in% c("Nonresponse", "Response")) %>% 
  dplyr::select(Group, Celltype_Major, Proportion) %>%
  tidyr::pivot_wider(names_from = Group, values_from = Proportion) %>%
  dplyr::mutate(Fold_Change = Nonresponse / Response,
                Log2FC = log2(Nonresponse / Response)) %>%
  ggplot2::ggplot(aes(x = Celltype_Major, y = Log2FC, fill = Log2FC > 0)) +
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
ggsave(plot = Log2FC_Nonresponse_Response, filename = "5_ScRNA_Celltype_Major/Log2FC_Nonresponse_Response.pdf", 
       units = "in", width = 10, height = 6, dpi = 300)

# Save ####
save(Seurat_Celltype_Major, file = "1_ScRNA_SeuratData/Seurat_Celltype_Major.RData")
# load(file = "1_ScRNA_SeuratData/Seurat_Celltype_Major.RData")

# File Tree ####
fs::dir_tree("5_ScRNA_Celltype_Major/", max_depth = 2)

# End ####
sessionInfo()

