# Set ####
options(timeout = 7355608); getOption("timeout")
options(stringsAsFactors = FALSE)
getwd(); rm(list = ls()); gc()
set.seed(19981108)
if (!dir.exists("6_ScRNA_hdWGCNA_T_Cell")) {dir.create("6_ScRNA_hdWGCNA_T_Cell")}

# Library ####
library(Seurat); packageVersion("Seurat")
library(hdWGCNA); packageVersion("hdWGCNA")
library(tidyverse); packageVersion("tidyverse")
library(future); packageVersion("future")
library(harmony); packageVersion("harmony")
library(patchwork); packageVersion("patchwork")
library(ggraph); packageVersion("ggraph")
library(cowplot); packageVersion("cowplot")
library(clustree); packageVersion("clustree")
library(ggrepel); packageVersion("ggrepel")
library(ggrastr); packageVersion("ggrastr")
library(BPCells); packageVersion("BPCells")
library(ggpubr); packageVersion("ggpubr")
library(ggradar); packageVersion("ggradar")
library(cols4all); packageVersion("cols4all")
library(enrichR); packageVersion("enrichR")
library(GeneOverlap); packageVersion("GeneOverlap")
library(fgsea); packageVersion("fgsea")

# Load ####
load(file = "1_ScRNA_SeuratData/Seurat_Celltype_T_Cell.RData")
Seurat_Celltype_T_Cell
head(Seurat_Celltype_T_Cell@meta.data, n = 10)
DefaultAssay(Seurat_Celltype_T_Cell) <- "RNA"
Seurat_hdWGCNA_T_Cell <- Seurat_Celltype_T_Cell; rm(Seurat_Celltype_T_Cell); gc()

# [Step 1] SetupForWGCNA ####
Seurat_hdWGCNA_T_Cell <- SetupForWGCNA(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  gene_select = "fraction", # fraction & variable
  fraction = 0.05,
  assay = "RNA", 
  wgcna_name = "T_Cell")

# [Step 2] MetacellsByGroups ####
Seurat_hdWGCNA_T_Cell <- MetacellsByGroups(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  group.by = c("Sample", "Celltype_T_Cell"), 
  reduction = "harmony", 
  k = 25, max_shared = 10,
  ident.group = "Celltype_T_Cell",
  wgcna_name = "T_Cell")
length(unique(Seurat_hdWGCNA_T_Cell$metacell_grouping))
head(Seurat_hdWGCNA_T_Cell@meta.data[, c("Sample", "Celltype_T_Cell", "metacell_grouping")], n = 10)
table(Seurat_hdWGCNA_T_Cell$metacell_grouping)

# [Step 3] NormalizeMetacells ####
Seurat_hdWGCNA_T_Cell <- NormalizeMetacells(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  wgcna_name = "T_Cell")
names(Seurat_hdWGCNA_T_Cell@misc)
names(Seurat_hdWGCNA_T_Cell@misc$T_Cell)
Seurat_hdWGCNA_T_Cell@misc$T_Cell$wgcna_metacell_obj
dim(GetAssayData(Seurat_hdWGCNA_T_Cell@misc$T_Cell$wgcna_metacell_obj,
                 assay = "RNA", 
                 layer = "data"))

# [Step 4] SetDatExpr ####
Seurat_hdWGCNA_T_Cell <- SetDatExpr(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  assay = "RNA", 
  slot = "data")

# [Step 5] TestSoftPowers ####
Seurat_hdWGCNA_T_Cell <- TestSoftPowers(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  networkType = "signed", 
  corFnc = "bicor", 
  wgcna_name = "T_Cell")
Plot_SoftPower <- PlotSoftPowers(seurat_obj = Seurat_hdWGCNA_T_Cell, wgcna_name = "T_Cell")
Plot_SoftPower <- wrap_plots(Plot_SoftPower, ncol = 2); Plot_SoftPower
ggsave(plot = Plot_SoftPower, filename = "6_ScRNA_hdWGCNA_T_Cell/Plot_SoftPower.pdf", 
       width = 6, height = 6, dpi = 300, units = "in")
SoftPower <- GetPowerTable(Seurat_hdWGCNA_T_Cell)
head(SoftPower, n = 10)
write_csv(SoftPower, file = "6_ScRNA_hdWGCNA_T_Cell/SoftPower.csv")

# [Step 6] ConstructNetwork ####
Seurat_hdWGCNA_T_Cell <- ConstructNetwork(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  wgcna_name = "T_Cell",
  tom_outdir = "6_ScRNA_hdWGCNA_T_Cell/TOM",
  tom_name = "T_Cell_TOM")
names(Seurat_hdWGCNA_T_Cell@misc$T_Cell)
head(Seurat_hdWGCNA_T_Cell@misc$T_Cell$wgcna_modules, n = 10)
pdf(file = "6_ScRNA_hdWGCNA_T_Cell/Plot_hdWGCNA_Dendrogram.pdf",
    width = 10,
    height = 6)
PlotDendrogram(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  wgcna_name = "T_Cell",
  main = "hdWGCNA Dendrogram")
dev.off()

# [Step 7] ModuleEigengenes ####
Seurat_hdWGCNA_T_Cell <- ModuleEigengenes(
  seurat_obj = Seurat_hdWGCNA_T_Cell, 
  group.by.vars = "Sample", 
  wgcna_name = "T_Cell")
ME <- GetMEs(seurat_obj = Seurat_hdWGCNA_T_Cell,
             wgcna_name = "T_Cell")
head(ME, n = 10)
write_csv(ME, file = "6_ScRNA_hdWGCNA_T_Cell/ME.csv")
save(ME, file = "6_ScRNA_hdWGCNA_T_Cell/ME.RData")

# [Step 8] ModuleConnectivity ####
Seurat_hdWGCNA_T_Cell <- ModuleConnectivity(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  wgcna_name = "T_Cell")
names(Seurat_hdWGCNA_T_Cell@misc$T_Cell$wgcna_modules)

# [Step 9] PlotKMEs ####
Plot_KME <- PlotKMEs(Seurat_hdWGCNA_T_Cell, ncol = 3, wgcna_name = "T_Cell"); Plot_KME
ggsave(plot = Plot_KME, filename = "6_ScRNA_hdWGCNA_T_Cell/Plot_KME.pdf", 
       width = 9, height = 9, dpi = 300, units = "in")

# [Step 10] GetModules ####
Module <- GetModules(
  seurat_obj = Seurat_hdWGCNA_T_Cell, 
  wgcna_name = "T_Cell")
head(Module, n = 10); table(Module$module)
write_csv(Module, file = "6_ScRNA_hdWGCNA_T_Cell/Module.csv")
save(Module, file = "6_ScRNA_hdWGCNA_T_Cell/Module.RData")

# [Step 11] GetHubGenes ####
HubGene <- GetHubGenes(
  seurat_obj = Seurat_hdWGCNA_T_Cell, 
  n_hubs = 10, 
  wgcna_name = "T_Cell")
head(HubGene, n = 10)
write_csv(HubGene, file = "6_ScRNA_hdWGCNA_T_Cell/HubGene.csv")
save(HubGene, file = "6_ScRNA_hdWGCNA_T_Cell/HubGene.RData")

# [Step 12] ModuleExprScore ####
Seurat_hdWGCNA_T_Cell <- ModuleExprScore(
  seurat_obj = Seurat_hdWGCNA_T_Cell, 
  n_genes = 25, 
  method = "UCell", # Seurat [AddModuleScore]
  wgcna_name = "T_Cell")

# [Step 13] ModuleFeaturePlot ####
Plot_ModuleFeature <- ModuleFeaturePlot(
  seurat_obj = Seurat_hdWGCNA_T_Cell, 
  features = "hMEs", # plot the hMEs
  order = TRUE, # order so the points with highest hMEs are on top
  reduction = "harmony_umap",
  wgcna_name = "T_Cell")
Plot_ModuleFeature <- wrap_plots(Plot_ModuleFeature, ncol = 3)
Plot_ModuleFeature
ggsave(plot = Plot_ModuleFeature, filename = "6_ScRNA_hdWGCNA_T_Cell/Plot_ModuleFeature.pdf", 
       width = 12, height = 12, dpi = 300, units = "in")

Plot_ModuleFeature_Score <- ModuleFeaturePlot(
  seurat_obj = Seurat_hdWGCNA_T_Cell, 
  features = "scores", # plot the hMEs
  order = "shuffle", # order so the points with highest hMEs are on top
  reduction = "harmony_umap",
  ucell = TRUE,
  wgcna_name = "T_Cell")
Plot_ModuleFeature_Score <- wrap_plots(Plot_ModuleFeature_Score, ncol = 3)
Plot_ModuleFeature_Score
ggsave(plot = Plot_ModuleFeature_Score, filename = "6_ScRNA_hdWGCNA_T_Cell/Plot_ModuleFeature_Score.pdf", 
       width = 12, height = 12, dpi = 300, units = "in")

# [Step 14] ModuleCorrelogram ####
pdf(file = "6_ScRNA_hdWGCNA_T_Cell/Plot_ModuleCorrelogram.pdf",
    width = 8,
    height = 8)
ModuleCorrelogram(
  seurat_obj = Seurat_hdWGCNA_T_Cell, 
  wgcna_name = "T_Cell")
dev.off()

# [Step 15] Dot_hdWGCNA ####
Seurat_hdWGCNA_T_Cell@meta.data <- cbind(Seurat_hdWGCNA_T_Cell@meta.data, ME)
head(Seurat_hdWGCNA_T_Cell@meta.data, n = 10)
Plot_Dot_hdWGCNA <- ggplot(
  data = DotPlot(
    object = Seurat_hdWGCNA_T_Cell,
    group.by = "Celltype_T_Cell",
    features = levels(Module$module))$data,
  mapping = aes(
    x = features.plot, y = id,
    size = pct.exp, fill = avg.exp.scaled)) +
  geom_point(shape = 21, color = "#000000FF") +
  scale_size_continuous(
    range = c(0, 8), name = "% Expressed",
    breaks = c(25, 50, 75, 100),
    limits = c(0, 100)) +
  scale_fill_gradient2(
    low = "#5DB1DDFF", high = "#CC0000FF",
    midpoint = 0, name = "Module Score") +
  facet_grid(
    . ~ features.plot,
    scales = "free_x", space = "free_x") +
  theme_classic() +
  labs(title = "hdWGCNA Module Activity in Celltype",
    x = "", y = "") +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    strip.text = element_text(size = 10, face = "bold", angle = 90, hjust = 0),
    strip.background = element_rect(fill = "#FFFFFFFF", color = "#FFFFFFFF"),
    legend.position = "right"); Plot_Dot_hdWGCNA
ggsave(plot = Plot_Dot_hdWGCNA, 
       filename = "6_ScRNA_hdWGCNA_T_Cell/Plot_Dot_hdWGCNA.pdf", 
       width = 6, height = 6, dpi = 300, units = "in")

# [Step 16] ModuleRadarPlot ####
Plot_ModuleRadar_Group <- ModuleRadarPlot(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  group.by = "Group",
  barcodes = Seurat_hdWGCNA_T_Cell@meta.data %>% rownames(),
  axis.label.size = 4,
  grid.label.size = 4); Plot_ModuleRadar_Group
ggsave(plot = Plot_ModuleRadar_Group, 
       filename = "6_ScRNA_hdWGCNA_T_Cell/Plot_ModuleRadar_Group.pdf", 
       width = 12, height = 12, dpi = 300, units = "in")

Plot_ModuleRadar_Celltype <- ModuleRadarPlot(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  group.by = "Celltype_T_Cell",
  barcodes = Seurat_hdWGCNA_T_Cell@meta.data %>% rownames(),
  axis.label.size = 4,
  grid.label.size = 4); Plot_ModuleRadar_Celltype
ggsave(plot = Plot_ModuleRadar_Celltype, 
       filename = "6_ScRNA_hdWGCNA_T_Cell/Plot_ModuleRadar_Celltype.pdf", 
       width = 12, height = 12, dpi = 300, units = "in")

# [Step 17] ModuleNetworkPlot ####
table(Seurat_hdWGCNA_T_Cell@misc$T_Cell$wgcna_modules$color)
ModuleNetworkPlot(
  seurat_obj = Seurat_hdWGCNA_T_Cell, 
  outdir = "6_ScRNA_hdWGCNA_T_Cell/ModuleNetworkPlot", 
  n_inner = 20,
  n_outer = 30,
  n_conns = Inf,
  n_genes = 30, 
  plot_size = c(10, 10),
  edge_threshold = 0.02,  
  vertex.label.cex = 1, 
  layout = "fr", 
  wgcna_name = "T_Cell")

# [Step 18] HubGeneNetworkPlot ####
pdf(file = "6_ScRNA_hdWGCNA_T_Cell/Plot_HubGeneNetwork.pdf",
    width = 5,
    height = 5)
options(future.globals.maxSize = 40 * 1024^3)
HubGeneNetworkPlot(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  n_hubs = 2, n_other = 5,
  edge_prop = 0.75,
  mods = "all",
  wgcna_name = "T_Cell")
dev.off()

head(Module, n = 10)
levels(Module$module)
HubGeneNetworkPlot(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  n_hubs = 2, n_other = 5,
  edge_prop = 0.75,
  mods = c("pink", "blue", "turquoise", "brown", 
           "green", "yellow", "black"),
  wgcna_name = "T_Cell")

# [Step 19] ModuleUMAPPlot ####
Seurat_hdWGCNA_T_Cell <- RunModuleUMAP(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors = 15, # neighbors parameter for UMAP
  min_dist = 0.1, # min distance between points in UMAP space
  wgcna_name = "T_Cell")
UMAP <- GetModuleUMAP(Seurat_hdWGCNA_T_Cell); head(UMAP, n = 10)
ggplot(UMAP, aes(x = UMAP1, y= UMAP2)) +
  geom_point(color = UMAP$color, size = UMAP$kME*2) +
  umap_theme()

pdf(file = "6_ScRNA_hdWGCNA_T_Cell/Plot_ModuleUMAP.pdf",
    width = 10,
    height = 10)
ModuleUMAPPlot(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  edge.alpha = 0.25,
  sample_edges = TRUE,
  edge_prop = 0.1, # proportion of edges to sample (20% here)
  label_hubs = 2 ,# how many hub genes to plot per module?
  keep_grey_edges = FALSE, 
  wgcna_name = "T_Cell")
dev.off()

# [Step 20] ModuleTraitCorrelation ####
Seurat_hdWGCNA_T_Cell@meta.data <- Seurat_hdWGCNA_T_Cell@meta.data %>% 
  dplyr::mutate(Trait = dplyr::case_when(
    Seurat_hdWGCNA_T_Cell@meta.data$Group == "Nonresponse" ~ "1", 
    Seurat_hdWGCNA_T_Cell@meta.data$Group == "Response" ~ "0", 
    TRUE ~ "Unknown"))
Seurat_hdWGCNA_T_Cell$Trait <- as.numeric(Seurat_hdWGCNA_T_Cell$Trait)
Seurat_hdWGCNA_T_Cell$Cell_Cycle <- as.factor(Seurat_hdWGCNA_T_Cell$Phase)
Seurat_hdWGCNA_T_Cell <- ModuleTraitCorrelation(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  traits = c("Trait", "Cell_Cycle"),
  group.by = "Celltype_T_Cell", 
  wgcna_name = "T_Cell")
Trait_Correlation <- GetModuleTraitCorrelation(Seurat_hdWGCNA_T_Cell)
names(Trait_Correlation)
names(Trait_Correlation$cor)

pdf(file = "6_ScRNA_hdWGCNA_T_Cell/Plot_ModuleTraitCorrelation.pdf",
    width = 8,
    height = 8)
PlotModuleTraitCorrelation(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  label = "fdr",
  label_symbol = "stars",
  text_size = 10,
  text_digits = 10,
  text_color = "#000000FF",
  high_color = "#CC0000FF",
  mid_color = "#FFFFFFFF",
  low_color = "#5DB1DDFF",
  plot_max = 1,
  combine = TRUE, 
  wgcna_name = "T_Cell")
dev.off()

# [Step 21] RunEnrichr ####
Pathway <- c(
  "GO_Biological_Process_2023", 
  "GO_Cellular_Component_2023", 
  "GO_Molecular_Function_2023", 
  "KEGG_2021_Human")
Seurat_hdWGCNA_T_Cell <- RunEnrichr(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  dbs = Pathway,
  max_genes = Inf, 
  wgcna_name = "T_Cell")
Enrich <- GetEnrichrTable(Seurat_hdWGCNA_T_Cell)
head(Enrich, n = 10)
write_csv(Enrich, file = "6_ScRNA_hdWGCNA_T_Cell/Enrich.csv")
save(Enrich, file = "6_ScRNA_hdWGCNA_T_Cell/Enrich.RData")

EnrichrBarPlot(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  outdir = "6_ScRNA_hdWGCNA_T_Cell/Plot_Enrichr", 
  p_cutoff = 0.05,
  p_adj = TRUE,
  n_terms = 10, # number of enriched terms to show (sometimes more are shown if there are ties)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale = TRUE, # do you want to show the enrichment as a log scale?
  wgcna_name = "T_Cell")

Plot_EnrichrDot <- EnrichrDotPlot(
  seurat_obj = Seurat_hdWGCNA_T_Cell,
  mods = "all", # use all modules (default)
  database = "KEGG_2021_Human", # this must match one of the dbs used previously
  n_terms = 10, # number of terms per module
  term_size = 8, # font size for the terms
  p_adj = TRUE, # show the p-val or adjusted p-val?
  wgcna_name = "T_Cell") + 
  scale_color_stepsn(colors=rev(viridis::magma(256))) + 
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    legend.position = "right"); Plot_EnrichrDot
ggsave(plot = Plot_EnrichrDot, 
       filename = "6_ScRNA_hdWGCNA_T_Cell/Plot_EnrichrDot.pdf", 
       width = 8, height = 10, dpi = 300, units = "in")

# [Step 22] Save ####
save(Seurat_hdWGCNA_T_Cell, file = "1_ScRNA_SeuratData/Seurat_hdWGCNA_T_Cell.RData")
# load(file = "1_ScRNA_SeuratData/Seurat_hdWGCNA_T_Cell.RData")

# File Tree ####
fs::dir_tree("6_ScRNA_hdWGCNA_T_Cell/", max_depth = 2)

# End ####
sessionInfo()
