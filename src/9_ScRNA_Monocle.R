# Set ####
options(timeout = 7355608); getOption("timeout")
options(stringsAsFactors = FALSE)
getwd(); rm(list = ls()); gc()
set.seed(19981108)
if (!dir.exists("9_ScRNA_Monocle")) {dir.create("9_ScRNA_Monocle")}

# Library ####
library(Seurat); packageVersion("Seurat")
library(tidyverse); packageVersion("tidyverse")
library(patchwork); packageVersion("patchwork")
library(ggrepel); packageVersion("ggrepel")
library(ggrastr); packageVersion("ggrastr")
library(ggpubr); packageVersion("ggpubr")
library(monocle); packageVersion("monocle")
library(igraph); packageVersion("igraph")
library(cowplot); packageVersion("cowplot")

# Load ####
load(file = "1_ScRNA_SeuratData/Seurat_Celltype_Major.RData")
DefaultAssay(Seurat_Celltype_Major) <- "RNA"; DefaultAssay(Seurat_Celltype_Major)
Seurat_Monocle <- subset(Seurat_Celltype_Major, subset = Celltype_Major == "T_Cell"); rm(Seurat_Celltype_Major); gc()
Idents(Seurat_Monocle) <- Seurat_Monocle@meta.data$Group
table(Idents(Seurat_Monocle))
table(Seurat_Monocle@meta.data$Celltype_T_Cell)

# Seurat => Monocle
Data_E <- as(as.matrix(Seurat_Monocle@assays$RNA@layers$counts), "sparseMatrix")
Data_P <- Seurat_Monocle@meta.data
Data_F <- data.frame(gene_short_name = rownames(Seurat_Monocle), 
                     row.names = rownames(Seurat_Monocle))

Data_P <- new("AnnotatedDataFrame", data = Data_P)
Data_F <- new("AnnotatedDataFrame", data = Data_F)
save(Data_E, file = "9_ScRNA_Monocle/Data_E.RData")
save(Data_P, file = "9_ScRNA_Monocle/Data_P.RData")
save(Data_F, file = "9_ScRNA_Monocle/Data_F.RData")

load(file = "9_ScRNA_Monocle/Data_E.RData")
load(file = "9_ScRNA_Monocle/Data_P.RData")
load(file = "9_ScRNA_Monocle/Data_F.RData")

CDS <- monocle::newCellDataSet(cellData = Data_E, 
                               phenoData = Data_P, 
                               featureData = Data_F, 
                               lowerDetectionLimit = 0.5, 
                               expressionFamily = negbinomial.size())
CDS <- estimateSizeFactors(CDS)
CDS <- estimateDispersions(CDS)
CDS <- detectGenes(CDS, min_expr = 0.1)
print(head(fData(CDS)))
CDS <- CDS[fData(CDS)$num_cells_expressed > 10, ]
save(CDS, file = "9_ScRNA_Monocle/CDS.RData")

# differentialGeneTest ####
DGT <- differentialGeneTest(CDS,
                            verbose = TRUE,
                            fullModelFormulaStr = "~ Group", 
                            cores = 16)
save(DGT, file = "9_ScRNA_Monocle/DGT.RData")
head(DGT)
Order_Gene <- DGT %>%
  dplyr::filter(qval < 0.01) %>%
  dplyr::arrange(qval) %>%
  dplyr::pull(gene_short_name) %>%
  utils::head(2000)
Order_Gene

# setOrderingFilter ####
CDS <- monocle::setOrderingFilter(CDS, Order_Gene)
Plot_Order_Gene <- plot_ordering_genes(CDS) + 
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold")); Plot_Order_Gene
ggsave(plot = Plot_Order_Gene, filename = "9_ScRNA_Monocle/Plot_Order_Gene.pdf", 
       width = 6, height = 6, dpi = 300, units = "in")
save(CDS, file = "9_ScRNA_Monocle/CDS.RData")

# reduceDimension & orderCells ####
CDS <- monocle::reduceDimension(cds = CDS, 
                                max_components = 2, 
                                verbose = TRUE,
                                method = "DDRTree")
save(CDS, file = "9_ScRNA_Monocle/CDS.RData")
load(file = "9_ScRNA_Monocle/CDS.RData")

CDS <- monocle::orderCells(CDS)
CDS <- monocle::orderCells(CDS, root_state = 7)
save(CDS, file = "9_ScRNA_Monocle/CDS.RData")
load(file = "9_ScRNA_Monocle/CDS.RData")

# Plot Pseudotime ####
# Pseudotime
plot_cell_trajectory(CDS, color_by = "State")
Plot_Pseudotime <- plot_cell_trajectory(CDS, color_by = "Pseudotime") + 
  theme_classic(base_size = 8) + 
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.key.size = unit(0.5, "in"), 
        legend.position = "right"); Plot_Pseudotime
ggsave(plot = Plot_Pseudotime, filename = "9_ScRNA_Monocle/Plot_Pseudotime.pdf", 
       width = 6, height = 5, dpi = 300, units = "in")

# Celltype
Plot_Pseudotime_Celltype <- plot_cell_trajectory(CDS, color_by = "Celltype_T_Cell") + 
  theme_classic(base_size = 8) + 
  facet_wrap("Celltype_T_Cell~ Group") + 
  guides(color = guide_legend(override.aes = list(size = 12))) +
  theme(panel.grid = element_blank(), 
        strip.background = element_rect(fill = "#FFFFFFFF", color = "#FFFFFFFF"), 
        strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.key.size = unit(0.5, "in"), 
        legend.position = "top"); Plot_Pseudotime_Celltype
ggsave(plot = Plot_Pseudotime_Celltype, filename = "9_ScRNA_Monocle/Plot_Pseudotime_Celltype.pdf", 
       width = 24, height = 20, dpi = 300, units = "in")

# State [Group] 
Plot_Pseudotime_State <- plot_cell_trajectory(CDS, color_by = "State") + 
  facet_wrap("~ Group", nrow = 1) + 
  theme_classic(base_size = 8) + 
  guides(color = guide_legend(override.aes = list(size = 10))) +
  # scale_color_manual(values = c("#FF1C1CFF", "#FB8B05FF", "#FFFF00FF", "#66FF66FF", 
  #                               "#1CCEFFFF", "#1CFFD2FF", "#A246FFFF", "#FF7AC2FF")) +
  scale_color_manual(values = c("#A246FFFF", "#1CFFD2FF", "#1CCEFFFF", 
                                "#66FF66FF", "#FFFF00FF", "#FB8B05FF", "#FF1C1CFF")) +
  theme(panel.grid = element_blank(), 
        strip.background = element_rect(fill = "#FFFFFFFF", color = "#FFFFFFFF"), 
        strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.key.size = unit(0.5, "in"), 
        legend.position = "right"); Plot_Pseudotime_State
ggsave(plot = Plot_Pseudotime_State, filename = "9_ScRNA_Monocle/Plot_Pseudotime_State.pdf", 
       width = 12, height = 5, dpi = 300, units = "in")

# Gene Biomarker
Plot_Gene_Pseudotime_Compare <- function(
  CDS, Gene, Group_Var = "Group",
  Group_Levels = c("Nonresponse", "Response"),
  Color_By = "State", Y_Lim = c(0, 120)) {
  Plot_List <- lapply(Group_Levels, function(Group_Name) {
    CDS_Gene <- CDS[Gene, ]
    CDS_Sub <- CDS_Gene[, pData(CDS_Gene)[[Group_Var]] == Group_Name]
    plot_genes_in_pseudotime(CDS_Sub, color_by = Color_By) +
      theme_classic(base_size = 8) +
      scale_y_continuous(limits = Y_Lim) +
      guides(color = guide_legend(override.aes = list(size = 10))) +
      scale_color_manual(values = c(
        "#FF1C1CFF", "#FB8B05FF", "#FFFF00FF", "#66FF66FF",
        "#1CCEFFFF", "#1CFFD2FF", "#A246FFFF", "#FF7AC2FF")) +
      theme(
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "#FFFFFFFF", color = "#FFFFFFFF"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.position = "none") +
      ggtitle(paste(Gene, Group_Name, sep = " - "))})
  Plot_List[[1]] | Plot_List[[2]]}
Plot_Pseudotime_Gene_ACTB <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "ACTB", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 130)); Plot_Pseudotime_Gene_ACTB
Plot_Pseudotime_Gene_ACTG1 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "ACTG1", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 50)); Plot_Pseudotime_Gene_ACTG1
Plot_Pseudotime_Gene_AKAP13 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "AKAP13", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 25)); Plot_Pseudotime_Gene_AKAP13
Plot_Pseudotime_Gene_CORO1A <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "CORO1A", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 25)); Plot_Pseudotime_Gene_CORO1A
Plot_Pseudotime_Gene_DPYD <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "DPYD", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 45)); Plot_Pseudotime_Gene_DPYD
Plot_Pseudotime_Gene_DUSP1 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "DUSP1", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 65)); Plot_Pseudotime_Gene_DUSP1
Plot_Pseudotime_Gene_FAM172A <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "FAM172A", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 40)); Plot_Pseudotime_Gene_FAM172A
Plot_Pseudotime_Gene_FOXN3 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "FOXN3", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 30)); Plot_Pseudotime_Gene_FOXN3
Plot_Pseudotime_Gene_IFITM1 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "IFITM1", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 60)); Plot_Pseudotime_Gene_IFITM1
Plot_Pseudotime_Gene_LRBA <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "LRBA", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 40)); Plot_Pseudotime_Gene_LRBA
Plot_Pseudotime_Gene_MYL12A <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "MYL12A", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 30)); Plot_Pseudotime_Gene_MYL12A
Plot_Pseudotime_Gene_PFN1 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "PFN1", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 100)); Plot_Pseudotime_Gene_PFN1
Plot_Pseudotime_Gene_PLAC8 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "PLAC8", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 40)); Plot_Pseudotime_Gene_PLAC8
Plot_Pseudotime_Gene_RNF213 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "RNF213", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 35)); Plot_Pseudotime_Gene_RNF213
Plot_Pseudotime_Gene_S100A4 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "S100A4", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 70)); Plot_Pseudotime_Gene_S100A4
Plot_Pseudotime_Gene_S100A6 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "S100A6", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 50)); Plot_Pseudotime_Gene_S100A6
Plot_Pseudotime_Gene_SERF2 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "SERF2", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 30)); Plot_Pseudotime_Gene_SERF2
Plot_Pseudotime_Gene_SH3KBP1 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "SH3KBP1", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 30)); Plot_Pseudotime_Gene_SH3KBP1
Plot_Pseudotime_Gene_SKAP1 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "SKAP1", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 50)); Plot_Pseudotime_Gene_SKAP1
Plot_Pseudotime_Gene_SP100 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "SP100", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 40)); Plot_Pseudotime_Gene_SP100
Plot_Pseudotime_Gene_TMSB4X <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "TMSB4X", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 180)); Plot_Pseudotime_Gene_TMSB4X
Plot_Pseudotime_Gene_TXNIP <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "TXNIP", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 80)); Plot_Pseudotime_Gene_TXNIP
Plot_Pseudotime_Gene_UTRN <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "UTRN", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 70)); Plot_Pseudotime_Gene_UTRN
Plot_Pseudotime_Gene_FOXP3 <-Plot_Gene_Pseudotime_Compare(
  CDS = CDS, Gene = "FOXP3", 
  Group_Levels = c("Nonresponse", "Response"), 
  Y_Lim = c(0, 10)); Plot_Pseudotime_Gene_FOXP3
Plot_Pseudotime_Gene_All <- cowplot::plot_grid(
  Plot_Pseudotime_Gene_ACTB, Plot_Pseudotime_Gene_ACTG1, Plot_Pseudotime_Gene_AKAP13, 
  Plot_Pseudotime_Gene_CORO1A, Plot_Pseudotime_Gene_DPYD, Plot_Pseudotime_Gene_DUSP1, 
  Plot_Pseudotime_Gene_FAM172A, Plot_Pseudotime_Gene_FOXN3, Plot_Pseudotime_Gene_IFITM1, 
  Plot_Pseudotime_Gene_LRBA, Plot_Pseudotime_Gene_MYL12A, Plot_Pseudotime_Gene_PFN1, 
  Plot_Pseudotime_Gene_PLAC8, Plot_Pseudotime_Gene_RNF213, Plot_Pseudotime_Gene_S100A4, 
  Plot_Pseudotime_Gene_S100A6, Plot_Pseudotime_Gene_SERF2, Plot_Pseudotime_Gene_SH3KBP1, 
  Plot_Pseudotime_Gene_SKAP1, Plot_Pseudotime_Gene_SP100, Plot_Pseudotime_Gene_TMSB4X, 
  Plot_Pseudotime_Gene_TXNIP, Plot_Pseudotime_Gene_UTRN, 
  ncol = 6); Plot_Pseudotime_Gene_All
ggsave(plot = Plot_Pseudotime_Gene_All, filename = "9_ScRNA_Monocle/Plot_Pseudotime_Gene_All.pdf", 
       width = 30, height = 20, dpi = 300, units = "in")

Plot_Pseudotime_Gene_Biomarker <- cowplot::plot_grid(
  Plot_Pseudotime_Gene_ACTB, Plot_Pseudotime_Gene_ACTG1, Plot_Pseudotime_Gene_AKAP13, 
  Plot_Pseudotime_Gene_CORO1A, Plot_Pseudotime_Gene_DUSP1, Plot_Pseudotime_Gene_MYL12A, 
  Plot_Pseudotime_Gene_PFN1, Plot_Pseudotime_Gene_S100A4, Plot_Pseudotime_Gene_S100A6,
  Plot_Pseudotime_Gene_SH3KBP1, Plot_Pseudotime_Gene_SKAP1, Plot_Pseudotime_Gene_TMSB4X, Plot_Pseudotime_Gene_UTRN, 
  ncol = 2); Plot_Pseudotime_Gene_Biomarker
Plot_Pseudotime_Gene_Biomarker <- cowplot::plot_grid(
  Plot_Pseudotime_Gene_ACTB, Plot_Pseudotime_Gene_ACTG1, 
  Plot_Pseudotime_Gene_PFN1, Plot_Pseudotime_Gene_S100A4, Plot_Pseudotime_Gene_S100A6,
  ncol = 1); Plot_Pseudotime_Gene_Biomarker
ggsave(plot = Plot_Pseudotime_Gene_Biomarker, filename = "9_ScRNA_Monocle/Plot_Pseudotime_Gene_Biomarker.pdf", 
       width = 8, height = 25, dpi = 300, units = "in", limitsize = FALSE)

# Time_Gene
Time_Gene <- differentialGeneTest(CDS, 
                                  fullModelFormulaStr = "~ sm.ns(Pseudotime)",
                                  cores = 16)
Time_Gene_Significant <- row.names(subset(Time_Gene, qval < 0.01))
Time_Gene_Significant
Pseudotime_Heatmap <- plot_pseudotime_heatmap(CDS[Time_Gene_Significant,],
                                              num_clusters = 5,
                                              cores = 16,
                                              show_rownames = F, 
                                              return_heatmap = T)
Pseudotime_Heatmap
ggsave(plot = Pseudotime_Heatmap, filename = "9_ScRNA_Monocle/Pseudotime_Heatmap.pdf", 
       width = 6, height = 18, dpi = 300, units = "in")
save(CDS, file = "9_ScRNA_Monocle/CDS.RData")

Pseudotime_Heatmap$tree_row
Clusters <- cutree(Pseudotime_Heatmap$tree_row, k = 5)
Clusters <- data.frame(Clusters)
Clusters[, 1] <- as.character(Clusters[, 1])
colnames(Clusters) <- "Gene_Clusters"
table(Clusters)
Clusters <- Clusters %>% 
  tibble::rownames_to_column(var = "Gene") %>% 
  tibble::as_tibble()
Clusters
save(Clusters, file = "9_ScRNA_Monocle/Clusters.RData")
load(file = "9_ScRNA_Monocle/GDT/Clusters.RData")

Cluster_Up <- Clusters %>% 
  dplyr::filter(Gene_Clusters %in% c("1", "2"))
Cluster_Up
Cluster_Down <- Clusters %>% 
  dplyr::filter(Gene_Clusters %in% c("3", "4", "5"))
Cluster_Down

# GO KEGG Up
Gene_Entrez <- bitr(geneID = Cluster_Up$Gene, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# GO 
Analyze_GO <- enrichGO(gene = Gene_Entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "ALL",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       pAdjustMethod = "fdr")

# KEGG 
Analyze_KEGG <- enrichKEGG(gene = Gene_Entrez$ENTREZID, 
                           organism = "hsa",
                           keyType = "kegg",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           pAdjustMethod = "fdr")

# Result
Result_GO <- Analyze_GO %>% 
  base::as.data.frame() %>% 
  tibble::as_tibble()
Result_KEGG <- Analyze_KEGG %>% 
  base::as.data.frame() %>% 
  tibble::as_tibble() %>% 
  dplyr::filter(!is.na(category))
Result_GO <- setReadable(Analyze_GO, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") %>% 
  base::as.data.frame() %>% 
  tibble::as_tibble()
Result_KEGG <- setReadable(Analyze_KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") %>% 
  base::as.data.frame() %>% 
  tibble::as_tibble() %>% 
  dplyr::filter(!is.na(category))

Result_GO <- Result_GO %>% 
  dplyr::mutate(ONTOLOGY = as.character(ONTOLOGY)) %>% 
  dplyr::group_by(ONTOLOGY) %>%
  dplyr::arrange(Count, .by_group = TRUE) %>%
  dplyr::ungroup()
Result_KEGG <- Result_KEGG %>%
  dplyr::arrange(Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG")
Result_All <- bind_rows(Result_GO %>% 
                          dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "CC", "MF"))), 
                        Result_KEGG %>%  
                          dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = "KEGG"))) %>%
  dplyr::arrange(ONTOLOGY, Count) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))
table(Result_All$ONTOLOGY)
write_csv(x = Result_All, file = "9_ScRNA_Monocle/GO_KEGG_Cluster_Up.csv")
Plot_Data <- Result_All %>% 
  dplyr::group_by(ONTOLOGY) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()
table(Plot_Data$ONTOLOGY)

# Plot
Plot_GO_KEGG <- ggplot() +
  # BP
  geom_point(data = Plot_Data %>% filter(ONTOLOGY == "BP"),
             mapping = aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), 
             shape = 21) + 
  scale_fill_gradient(low = "#8C96C6", high = "#8C6BB1", name = "BP p.adjust") + 
  ggnewscale::new_scale_fill() + 
  # CC
  geom_point(data = Plot_Data %>% filter(ONTOLOGY == "CC"), 
             mapping = aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), 
             shape = 21) + 
  scale_fill_gradient(low = "#FDD49E", high = "#D7301F", name = "CC p.adjust") + 
  ggnewscale::new_scale_fill() +
  # MF 
  geom_point(data = Plot_Data %>% filter(ONTOLOGY == "MF"),
             mapping = aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), 
             shape = 21) + 
  scale_fill_gradient(low = "#A6BDDB", high = "#0570B0", name = "MF p.adjust") + 
  ggnewscale::new_scale_fill() +
  # KEGG
  geom_point(data = Plot_Data %>% filter(ONTOLOGY == "KEGG"), 
             mapping = aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), 
             shape = 21) +
  scale_fill_gradient(low = "#A1D99B", high = "#238B45", name = "KEGG p.adjust") + 
  ggnewscale::new_scale_fill() +
  guides(y = "axis_nested", 
         y.sec = guide_axis_manual(breaks = seq_len(nrow(Plot_Data)),  
                                   labels = Plot_Data$Description), 
         fill1 = guide_legend(order = 1, title = "KEGG p.adjust"),
         fill2 = guide_legend(order = 2, title = "MF p.adjust"),
         fill3 = guide_legend(order = 3, title = "CC p.adjust"),
         fill4 = guide_legend(order = 4, title = "BP p.adjust"),
         size = guide_legend(order = 5, title = "Gene Count")) +
  ggtitle(label = "GO and KEGG Pathway Enrichment") + 
  labs(x = "Gene Count", y = "Description") + 
  scale_size(range = c(2, 10),
             guide = guide_legend(override.aes = list(fill = "#000000"))) + 
  theme_bw() + 
  theme(ggh4x.axis.nestline.y = element_line(linewidth = 5, color = c("#88419D", "#F46D43", "#41B6C4", "#74C476")),
        ggh4x.axis.nesttext.y = element_text(colour = c("#88419D", "#F46D43", "#41B6C4", "#74C476")),
        legend.position = "left", 
        legend.background = element_rect(color = "#969696", linewidth = 0.5),
        panel.border = element_rect(linewidth = 0.5),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
        axis.text = element_text(color = "#000000", size = 11),
        axis.text.y = element_text(color = rep(c("#88419D", "#FC4E2A", "#41B6C4", "#74C476"), times = c(10, 10, 10, 10))),
        axis.text.y.left = element_blank(),
        axis.ticks.length.y.left = unit(10, "pt"),
        axis.ticks.y.left = element_line(color = NA),
        axis.title = element_text(color = "#000000", size = 15),
        plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)) + 
  coord_cartesian(clip = "off")
Plot_GO_KEGG
ggsave(filename = "9_ScRNA_Monocle/Plot_GO_KEGG_Cluster_Up.pdf", plot = Plot_GO_KEGG, 
       width = 12, height = 12, dpi = 300, units = "in")

# GO KEGG Down
Gene_Entrez <- bitr(geneID = Cluster_Down$Gene, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# GO 
Analyze_GO <- enrichGO(gene = Gene_Entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "ALL",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       pAdjustMethod = "fdr")

# KEGG 
Analyze_KEGG <- enrichKEGG(gene = Gene_Entrez$ENTREZID, 
                           organism = "hsa",
                           keyType = "kegg",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           pAdjustMethod = "fdr")

# Result
Result_GO <- Analyze_GO %>% 
  base::as.data.frame() %>% 
  tibble::as_tibble()
Result_KEGG <- Analyze_KEGG %>% 
  base::as.data.frame() %>% 
  tibble::as_tibble() %>% 
  dplyr::filter(!is.na(category))
Result_GO <- setReadable(Analyze_GO, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") %>% 
  base::as.data.frame() %>% 
  tibble::as_tibble()
Result_KEGG <- setReadable(Analyze_KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") %>% 
  base::as.data.frame() %>% 
  tibble::as_tibble() %>% 
  dplyr::filter(!is.na(category))

Result_GO <- Result_GO %>% 
  dplyr::mutate(ONTOLOGY = as.character(ONTOLOGY)) %>% 
  dplyr::group_by(ONTOLOGY) %>%
  dplyr::arrange(Count, .by_group = TRUE) %>%
  dplyr::ungroup()
Result_KEGG <- Result_KEGG %>%
  dplyr::arrange(Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG")
Result_All <- bind_rows(Result_GO %>% 
                          dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "CC", "MF"))), 
                        Result_KEGG %>%  
                          dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = "KEGG"))) %>%
  dplyr::arrange(ONTOLOGY, Count) %>%
  dplyr::mutate(Description = factor(Description, levels = unique(Description)))
table(Result_All$ONTOLOGY)
write_csv(x = Result_All, file = "9_ScRNA_Monocle/GO_KEGG_Cluster_Down.csv")
Plot_Data <- Result_All %>% 
  dplyr::group_by(ONTOLOGY) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()
table(Plot_Data$ONTOLOGY)

# Plot
Plot_GO_KEGG <- ggplot() +
  # BP
  geom_point(data = Plot_Data %>% filter(ONTOLOGY == "BP"),
             mapping = aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), 
             shape = 21) + 
  scale_fill_gradient(low = "#8C96C6", high = "#8C6BB1", name = "BP p.adjust") + 
  ggnewscale::new_scale_fill() + 
  # CC
  geom_point(data = Plot_Data %>% filter(ONTOLOGY == "CC"), 
             mapping = aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), 
             shape = 21) + 
  scale_fill_gradient(low = "#FDD49E", high = "#D7301F", name = "CC p.adjust") + 
  ggnewscale::new_scale_fill() +
  # MF 
  geom_point(data = Plot_Data %>% filter(ONTOLOGY == "MF"),
             mapping = aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), 
             shape = 21) + 
  scale_fill_gradient(low = "#A6BDDB", high = "#0570B0", name = "MF p.adjust") + 
  ggnewscale::new_scale_fill() +
  # KEGG
  geom_point(data = Plot_Data %>% filter(ONTOLOGY == "KEGG"), 
             mapping = aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), 
             shape = 21) +
  scale_fill_gradient(low = "#A1D99B", high = "#238B45", name = "KEGG p.adjust") + 
  ggnewscale::new_scale_fill() +
  guides(y = "axis_nested", 
         y.sec = guide_axis_manual(breaks = seq_len(nrow(Plot_Data)),  
                                   labels = Plot_Data$Description), 
         fill1 = guide_legend(order = 1, title = "KEGG p.adjust"),
         fill2 = guide_legend(order = 2, title = "MF p.adjust"),
         fill3 = guide_legend(order = 3, title = "CC p.adjust"),
         fill4 = guide_legend(order = 4, title = "BP p.adjust"),
         size = guide_legend(order = 5, title = "Gene Count")) +
  ggtitle(label = "GO and KEGG Pathway Enrichment") + 
  labs(x = "Gene Count", y = "Description") + 
  scale_size(range = c(2, 10),
             guide = guide_legend(override.aes = list(fill = "#000000"))) + 
  theme_bw() + 
  theme(ggh4x.axis.nestline.y = element_line(linewidth = 5, color = c("#88419D", "#F46D43", "#41B6C4", "#74C476")),
        ggh4x.axis.nesttext.y = element_text(colour = c("#88419D", "#F46D43", "#41B6C4", "#74C476")),
        legend.position = "left", 
        legend.background = element_rect(color = "#969696", linewidth = 0.5),
        panel.border = element_rect(linewidth = 0.5),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
        axis.text = element_text(color = "#000000", size = 11),
        axis.text.y = element_text(color = rep(c("#88419D", "#FC4E2A", "#41B6C4", "#74C476"), times = c(10, 10, 10, 10))),
        axis.text.y.left = element_blank(),
        axis.ticks.length.y.left = unit(10, "pt"),
        axis.ticks.y.left = element_line(color = NA),
        axis.title = element_text(color = "#000000", size = 15),
        plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)) + 
  coord_cartesian(clip = "off")
Plot_GO_KEGG
ggsave(filename = "9_ScRNA_Monocle/Plot_GO_KEGG_Cluster_Down.pdf", plot = Plot_GO_KEGG, 
       width = 12, height = 12, dpi = 300, units = "in")

# File Tree ####
fs::dir_tree("9_ScRNA_Monocle/", max_depth = 2)

# End ####
sessionInfo()