# Set ####
options(timeout = 7355608); getOption("timeout")
options(stringsAsFactors = FALSE)
getwd(); rm(list = ls()); gc()
set.seed(19981108)
if (!dir.exists("7_ScRNA_DEA")) {dir.create("7_ScRNA_DEA")}

# Library ####
library(Seurat); packageVersion("Seurat")
library(tidyverse); packageVersion("tidyverse")
library(patchwork); packageVersion("patchwork")
library(ggrepel); packageVersion("ggrepel")
library(ggrastr); packageVersion("ggrastr")
library(ggalluvial); packageVersion("ggalluvial")
library(ggpubr); packageVersion("ggpubr")
library(cols4all); packageVersion("cols4all")
library(ggVennDiagram); packageVersion("ggVennDiagram")
library(ggforce); packageVersion("ggforce")
library(clusterProfiler); packageVersion("clusterProfiler")
library(org.Hs.eg.db); packageVersion("org.Hs.eg.db")
library(enrichplot); packageVersion("enrichplot")
library(ggh4x); packageVersion("ggh4x")
library(ggfun); packageVersion("ggfun")
library(ggnewscale); packageVersion("ggnewscale")
library(grid); packageVersion("grid")
library(gghalves); packageVersion("gghalves")
library(ggsankey); packageVersion("ggsankey")

# Load ####
load(file = "1_ScRNA_SeuratData/Seurat_Celltype_Major.RData")
head(Seurat_Celltype_Major@meta.data, n = 10)
table(Seurat_Celltype_Major@meta.data$Celltype_Major)
table(Seurat_Celltype_Major@meta.data$Group)
Idents(Seurat_Celltype_Major) <- "Group"
DefaultAssay(Seurat_Celltype_Major) <- "RNA"

# DEA ####
Celltype_DEA <- function(Celltype, Test_Method, Threshold_Log2FC, Threshold_PVal, Threshold_AdjPVal) {
  FindMarkers(
    object = subset(Seurat_Celltype_Major, 
                    subset = Celltype_Major == Celltype), 
    ident.1 = "Nonresponse", ident.2 = "Response", verbose = TRUE, only.pos = FALSE, 
    min.pct = 0.1, logfc.threshold = 0, test.use = Test_Method) %>% 
    tibble::rownames_to_column(var = "Gene_Symbol") %>% 
    dplyr::mutate(Celltype = Celltype) %>% 
    dplyr::filter(pct.1 > 0.2 & pct.2 > 0.2) %>% 
    dplyr::mutate(Regulation = case_when(avg_log2FC > Threshold_Log2FC ~ "Up", avg_log2FC < -Threshold_Log2FC ~ "Down", TRUE ~ "Not"),
                  Significant = case_when(p_val < Threshold_PVal ~ "Yes", TRUE ~ "No"), 
                  Significant_Adjusted = case_when(p_val_adj < Threshold_AdjPVal ~ "Yes", TRUE ~ "No"))}
DEG_T_Cell <- Celltype_DEA(
  Celltype = "T_Cell", 
  Test_Method = "wilcox", 
  Threshold_Log2FC = 0.6, 
  Threshold_PVal = 0.05, 
  Threshold_AdjPVal = 0.05) %>% 
  dplyr::rename(PValue = p_val, 
                Log2FoldChange = avg_log2FC, 
                PAdj = p_val_adj, 
                PCT1 = pct.1, 
                PCT2 = pct.2) %>% 
  dplyr::mutate(PAdj = ifelse(PAdj == 0, 1e-300, PAdj), 
                PAdj_Log = -log10(PAdj))
head(DEG_T_Cell, n = 10)
DEG_T_Cell_Sig <- DEG_T_Cell %>% 
  dplyr::filter(Significant == "Yes", 
                Significant_Adjusted == "Yes")
table(DEG_T_Cell_Sig$Regulation)
save(DEG_T_Cell, file = "7_ScRNA_DEA/DEG_T_Cell.RData")
write_csv(x = DEG_T_Cell, file = "7_ScRNA_DEA/DEG_T_Cell.csv")


# Plot_Volcano
Threshold_AdjPVal <- 0.05
Threshold_Log2FC <- 0.6
Colors_Volcano <-  c("#313695FF", "#4575B4FF", "#74ADD1FF", "#FDAE61FF", "#D73027FF", "#A50026FF") 

Plot_Volcano <- ggplot() +
  geom_point(data = DEG_T_Cell,
             mapping = aes(x = Log2FoldChange,
                           y = PAdj_Log,
                           colour = Log2FoldChange,
                           size = PAdj_Log)) + 
  geom_point(data = DEG_T_Cell %>% 
               dplyr::filter(Significant_Adjusted == "Yes") %>% 
               dplyr::arrange(desc(Log2FoldChange)) %>% 
               dplyr::slice(1:20),
             mapping = aes(x = Log2FoldChange,
                           y = PAdj_Log,
                           colour = Log2FoldChange,
                           size = PAdj_Log)) +
  geom_text_repel(data = DEG_T_Cell %>% 
                    dplyr::filter(Significant_Adjusted == "Yes") %>% 
                    dplyr::arrange(desc(Log2FoldChange)) %>% 
                    dplyr::slice(1:15) %>% 
                    dplyr::filter(Regulation == "Up"),
                  mapping = aes(x = Log2FoldChange,
                                y = PAdj_Log,
                                label = Gene_Symbol),
                  box.padding = 0.5,
                  nudge_x = 2,
                  nudge_y = 2,
                  segment.curvature = -0.1,
                  direction = "y",
                  force = 2,
                  max.overlaps = 30,
                  colour = "#F1441DFF") +
  geom_text_repel(data = DEG_T_Cell %>% 
                    dplyr::filter(Significant_Adjusted == "Yes") %>% 
                    dplyr::arrange(Log2FoldChange) %>% 
                    dplyr::slice(1:15) %>% 
                    dplyr::filter(Regulation == "Down"),
                  mapping = aes(x = Log2FoldChange,
                                y = PAdj_Log,
                                label = Gene_Symbol),
                  box.padding = 0.5,
                  nudge_x = -2,
                  nudge_y = 2,
                  segment.curvature = -0.1,
                  direction = "y",
                  force = 2,
                  max.overlaps = 30,
                  colour = "#2474B5FF") +
  scale_color_gradientn(name = "log2Fold Change",
                        colours = Colors_Volcano,
                        values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(name = "log2Fold Change",
                       colours = Colors_Volcano,
                       values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(-Threshold_Log2FC, Threshold_Log2FC), 
             linetype = 2, 
             linewidth = 1,
             colour = "#685E48FF") +
  geom_hline(yintercept = -log10(Threshold_AdjPVal), 
             linetype = 2, 
             linewidth = 1,
             colour = "#685E48FF") +
  scale_size(name = "- log10 Adjusted P-value", 
             range = c(1, 5)) +
  labs(title = "Volcano Plot of Differential Gene Expression", 
       subtitle = "Sensitive vs Resistant",
       x = "Log2 Fold Change", 
       y = "-log10 Adjusted P-value") + 
  xlim(c(-c(max(abs(DEG_T_Cell$Log2FoldChange))) - 2, 
         max(abs(DEG_T_Cell$Log2FoldChange)) + 2)) +  
  ylim(c(-20, max(abs(DEG_T_Cell$PAdj_Log)) + 100)) + 
  theme_bw()+
  theme(panel.grid = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        plot.subtitle = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.key.size = unit(0.5, "in")) +
  scale_size_continuous(range = c(0, 8)) +
  annotation_custom(grob = grid::segmentsGrob(y0 = unit(-10, "pt"),
                                              y1 = unit(-10, "pt"), 
                                              arrow = arrow(angle = 45, length = unit(.3, "cm"), ends = "first"),  
                                              gp = grid::gpar(lwd = 3, col = "#2474B5")), 
                    xmin = -c(max(abs(DEG_T_Cell$Log2FoldChange)) + 1), 
                    xmax = -1,
                    ymin = max(abs(DEG_T_Cell$PAdj_Log)) + 100, 
                    ymax = max(abs(DEG_T_Cell$PAdj_Log)) + 100) + 
  annotation_custom(grob = grid::textGrob(label = "Down Regulated Genes (68)", 
                                          gp = grid::gpar(col = "#2474B5")),
                    xmin = -c(max(abs(DEG_T_Cell$Log2FoldChange))), 
                    xmax = -2,
                    ymin = max(abs(DEG_T_Cell$PAdj_Log)) + 105, 
                    ymax = max(abs(DEG_T_Cell$PAdj_Log)) + 105) +
  
  annotation_custom(grob = grid::segmentsGrob(y0 = unit(-10, "pt"), 
                                              y1 = unit(-10, "pt"), 
                                              arrow = arrow(angle = 45, length = unit(.3, "cm"), ends = "last"), 
                                              gp = grid::gpar(lwd = 3, col = "#F1441D")), 
                    xmin = max(abs(DEG_T_Cell$Log2FoldChange) + 1), 
                    xmax = 1, 
                    ymin = max(abs(DEG_T_Cell$PAdj_Log)) + 100, 
                    ymax = max(abs(DEG_T_Cell$PAdj_Log)) + 100) + 
  annotation_custom(grob = grid::textGrob(label = "Up Regulated Genes (168)", 
                                          gp = grid::gpar(col = "#F1441D")), 
                    xmin = c(max(abs(DEG_T_Cell$Log2FoldChange))), 
                    xmax = 2, 
                    ymin = max(abs(DEG_T_Cell$PAdj_Log)) + 105, 
                    ymax = max(abs(DEG_T_Cell$PAdj_Log)) + 105); Plot_Volcano
ggsave(plot = Plot_Volcano, filename = "7_ScRNA_DEA/Plot_Volcano.pdf", 
       units = "in", width = 10, height = 10, dpi = 300)

# Module
load(file = "6_ScRNA_hdWGCNA_T_Cell/Module.RData")
head(Module, n = 10); table(Module$module)
Module_Select <- Module %>% 
  dplyr::filter(module %in% c("turquoise", "pink", "blue")) %>% 
  dplyr::select(c("gene_name", "module")) %>% 
  dplyr::rename(Gene_Symbol = gene_name, 
                Module = module)
head(Module_Select, n = 10); table(Module_Select$Module)

DEG_T_Cell_Module <- Module_Select %>% 
  dplyr::left_join(DEG_T_Cell, by = "Gene_Symbol") %>% 
  dplyr::filter(!is.na(Celltype)) %>% 
  dplyr::filter(Significant == "Yes") %>% 
  dplyr::filter(Significant_Adjusted == "Yes") %>% 
  dplyr::filter(PCT1 > 0.5) %>% 
  dplyr::filter(PCT2 > 0.5) %>% 
  dplyr::filter(abs(Log2FoldChange) > 0.6) %>% 
  dplyr::arrange(dplyr::desc(Log2FoldChange))
DEG_T_Cell_Module
write_csv(DEG_T_Cell_Module, file = "7_ScRNA_DEA/DEG_T_Cell_Module.csv")

# Venn ####
Venn_Data <- list(DEGs = DEG_T_Cell_Sig %>% 
                    dplyr::filter(!is.na(Celltype)) %>% 
                    dplyr::filter(Significant == "Yes") %>% 
                    dplyr::filter(Significant_Adjusted == "Yes") %>% 
                    dplyr::filter(PCT1 > 0.5) %>% 
                    dplyr::filter(PCT2 > 0.5) %>% 
                    dplyr::filter(abs(Log2FoldChange) > 0.6) %>% 
                    dplyr::pull(Gene_Symbol), 
                  Module_turquoise = Module_Select %>% 
                    dplyr::filter(Module == "turquoise") %>% 
                    dplyr::pull(Gene_Symbol), 
                  Module_pink = Module_Select %>% 
                    dplyr::filter(Module == "pink") %>% 
                    dplyr::pull(Gene_Symbol),  
                  Module_blue = Module_Select %>% 
                    dplyr::filter(Module == "blue") %>% 
                    dplyr::pull(Gene_Symbol))
Venn_Data <- list(DEGs = DEG_T_Cell_Sig %>% 
                    dplyr::filter(!is.na(Celltype)) %>% 
                    dplyr::filter(Significant == "Yes") %>% 
                    dplyr::filter(Significant_Adjusted == "Yes") %>% 
                    dplyr::filter(PCT1 > 0.5) %>% 
                    dplyr::filter(PCT2 > 0.5) %>% 
                    dplyr::filter(abs(Log2FoldChange) > 0.6) %>% 
                    dplyr::pull(Gene_Symbol), 
                  Module = Module_Select$Gene_Symbol)
Plot_Venn <- ggVennDiagram(Venn_Data, label = "count", edge_color = "#FFFFFFFF", label_color = "#FFFFFFFF") +  
  scale_fill_gradientn(colours = c("#F9D367FF", "#FF6961FF")) + 
  theme_void() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(title = "Venn Diagram of DEGs and Module Genes", fill = "Intersection\nCounts") + 
  coord_flip(); Plot_Venn
ggsave(filename = "7_ScRNA_DEA/Plot_Venn.pdf", plot = Plot_Venn, 
       width = 6, height = 4, dpi = 300, units = "in")

# Biomarker ####
MetaData_Celltype_Major <- Seurat_Celltype_Major@meta.data
MetaData_Celltype_Major[, DEG_T_Cell_Module$Gene_Symbol] <- FetchData(
  Seurat_Celltype_Major, assay = "RNA", vars = DEG_T_Cell_Module$Gene_Symbol)
head(MetaData_Celltype_Major, n = 10)
Plot_Box <- MetaData_Celltype_Major %>% 
  tibble::rownames_to_column(var = "Barcode") %>% 
  tibble::as_tibble() %>% 
  dplyr::filter(Celltype_Major == "T_Cell") %>% 
  dplyr::select(Barcode, all_of(DEG_T_Cell_Module$Gene_Symbol), Celltype_Major, Group) %>% 
  tidyr::pivot_longer(cols = -c(Barcode, Celltype_Major, Group),
                      names_to = "Gene", 
                      values_to = "Expression") %>% 
  # dplyr::filter(Expression != 0) %>% 
  ggplot2::ggplot(aes(x = Group, y = Expression, fill = Group)) +
  gghalves::geom_half_violin(
    side = "r", trim = FALSE,
    scale = "width", linewidth = 1, alpha = 0.8) +
  gghalves::geom_half_boxplot(
    side = "l", width = 0.6, outlier.shape = NA, linewidth = 1) +
  ggpubr::stat_compare_means(
    comparisons = list(c("Response", "Nonresponse")),
    method = "wilcox.test", label = "p.signif", size = 8, bracket.size = 1) +
  ggplot2::facet_wrap(~ Gene, nrow = 3) +
  ggplot2::coord_cartesian(ylim = c(0, 8)) +
  ggplot2::labs(title = "T Cell", y = "Normalized expression", x = "Group") +
  ggplot2::theme_classic(base_size = 10) + 
  ggplot2::theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"), 
    axis.title.y = element_text(size = 12, face = "bold"), 
    axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust = 1, vjust = 0.5), 
    axis.text.y = element_text(size = 12, face = "bold"), 
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"), 
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
    legend.key.size = unit(0.5, "in")) +
  ggplot2::scale_fill_manual(values = c("#FF6354FF", "#33CC33FF")); Plot_Box
ggsave(filename = "7_ScRNA_DEA/Plot_Box.pdf", plot = Plot_Box, 
       width = 10, height = 8, dpi = 300, units = "in")

Plot_Box <- MetaData_Celltype_Major %>% 
  tibble::rownames_to_column(var = "Barcode") %>% 
  tibble::as_tibble() %>% 
  dplyr::filter(Celltype_Major == "T_Cell") %>% 
  dplyr::select(Barcode, all_of(DEG_T_Cell_Module$Gene_Symbol), Celltype_T_Cell, Group) %>% 
  tidyr::pivot_longer(cols = -c(Barcode, Celltype_T_Cell, Group),
                      names_to = "Gene", 
                      values_to = "Expression") %>% 
  # dplyr::filter(Expression != 0) %>% 
  ggplot2::ggplot(aes(x = Group, y = Expression, fill = Group)) +
  gghalves::geom_half_violin(
    side = "r", trim = FALSE,
    scale = "width", linewidth = 1, alpha = 0.8) +
  gghalves::geom_half_boxplot(
    side = "l", width = 0.6, outlier.shape = NA, linewidth = 1) +
  ggpubr::stat_compare_means(
    comparisons = list(c("Response", "Nonresponse")),
    method = "wilcox.test", label = "p.signif", size = 8, bracket.size = 1) +
  ggplot2::facet_wrap(Celltype_T_Cell ~ Gene) +
  ggplot2::coord_cartesian(ylim = c(0, 8)) +
  ggplot2::labs(title = "T Cell", y = "Normalized expression", x = "Group") +
  ggplot2::theme_classic(base_size = 10) + 
  ggplot2::theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"), 
    axis.title.y = element_text(size = 12, face = "bold"), 
    axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust = 1, vjust = 0.5), 
    axis.text.y = element_text(size = 12, face = "bold"), 
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"), 
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
    legend.key.size = unit(0.5, "in")) +
  ggplot2::scale_fill_manual(values = c("#FF6354FF", "#33CC33FF")); Plot_Box
ggsave(filename = "7_ScRNA_DEA/Plot_Box_T_Cell.pdf", plot = Plot_Box, 
       width = 20, height = 30, dpi = 300, units = "in")

# GO KEGG
Gene_Entrez <- bitr(
  geneID = DEG_T_Cell_Module$Gene_Symbol, 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db)
Gene_Entrez_Universe <- bitr(
  geneID = DEG_T_Cell$Gene_Symbol, 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db)
Gene_Entrez_Universe <- unique(as.character(Gene_Entrez_Universe$ENTREZID))

# GO 
Analyze_GO <- enrichGO(
  gene = Gene_Entrez$ENTREZID, 
  OrgDb = org.Hs.eg.db, 
  keyType = "ENTREZID",
  ont = "ALL", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.2, 
  # universe = Gene_Entrez_Universe,
  pAdjustMethod = "fdr")

# KEGG 
Analyze_KEGG <- enrichKEGG(
  gene = Gene_Entrez$ENTREZID, 
  organism = "hsa", 
  keyType = "kegg", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.2, 
  # universe = Gene_Entrez_Universe,
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
write_csv(x = Result_All, file = "7_ScRNA_DEA/GO_KEGG_Result_All.csv")
Plot_Data <- Result_All %>% 
  dplyr::filter(
    Description %in% c(
      # BP
      "postsynaptic actin cytoskeleton organization",
      "regulation of synaptic vesicle endocytosis", 
      "regulation of actin cytoskeleton organization", 
      "actin filament organization", 
      "regulation of actomyosin structure organization", 
      "homotypic cell-cell adhesion", 
      "cell-substrate adhesion", 
      "chemotaxis", 
      # CC
      "postsynaptic actin cytoskeleton",
      "immunological synapse", 
      "cortical actin cytoskeleton", 
      "cortical cytoskeleton", 
      "actin filament", 
      "cell cortex",
      "cell leading edge", 
      "focal adhesion", 
      "cell-substrate junction",
      # MF
      "structural constituent of postsynaptic actin cytoskeleton", 
      "myosin heavy chain binding", 
      "structural constituent of postsynapse", 
      "structural constituent of synapse", 
      "myosin binding", 
      "actin monomer binding", 
      "actin binding", 
      # KEGG
      "Apoptosis", 
      "Adherens junction", 
      "Leukocyte transendothelial migration", 
      "Tight junction", 
      "Motor proteins", 
      "Focal adhesion", 
      "IgSF CAM signaling", 
      "Rap1 signaling pathway", 
      "Regulation of actin cytoskeleton"))
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
  scale_size(range = c(2, 10), guide = guide_legend(override.aes = list(fill = "#000000"))) + 
  theme_bw() + 
  theme(ggh4x.axis.nestline.y = element_line(linewidth = 5, color = c("#88419D", "#F46D43", "#41B6C4", "#74C476")),
        ggh4x.axis.nesttext.y = element_text(colour = c("#88419D", "#F46D43", "#41B6C4", "#74C476")),
        legend.position = "left", 
        legend.background = element_rect(color = "#969696", linewidth = 0.5),
        panel.border = element_rect(linewidth = 0.5),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
        axis.text = element_text(color = "#000000", size = 11),
        axis.text.y = element_text(color = rep(c("#88419D", "#FC4E2A", "#41B6C4", "#74C476"), times = c(8, 9, 7, 9))),
        axis.text.y.left = element_blank(),
        axis.ticks.length.y.left = unit(10, "pt"),
        axis.ticks.y.left = element_line(color = NA),
        axis.title = element_text(color = "#000000", size = 15),
        plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)) + 
  coord_cartesian(clip = "off")
Plot_GO_KEGG
ggsave(filename = "7_ScRNA_DEA/Plot_GO_KEGG.pdf", plot = Plot_GO_KEGG, 
       width = 10, height = 10, dpi = 300, units = "in")

Data_Sankey <- Plot_Data %>% 
  # dplyr::filter(ONTOLOGY == "KEGG") %>% 
  dplyr::select(Description, geneID, Count, ONTOLOGY) %>% 
  dplyr::mutate(geneID = str_split(geneID, "/")) %>%
  tidyr::unnest(geneID) %>%
  dplyr::rename(Gene_Symbol = geneID,
                Pathway = Description) %>% 
  ggsankey::make_long(Pathway, ONTOLOGY, Gene_Symbol) %>% 
  dplyr::mutate(node = factor(node))
Data_Sankey
Plot_Sankey <- ggplot(
  data = Data_Sankey, 
  mapping = aes(
    x = x, next_x = next_x, 
    node = node, next_node = next_node, 
    fill = node, label = node)) + 
  labs(title = "Receptor Ligands of Source (T Cell)") +
  geom_sankey(flow.alpha = 0.4, width = 0.12, smooth = 10) +
  geom_sankey_text(size = 4, type = "sankey", hjust = 0) +
  theme_void() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold")); Plot_Sankey
ggsave(plot = Plot_Sankey, filename = "7_ScRNA_DEA/Plot_Sankey.pdf",
       height = 8, width = 10, dpi = 300, units = "in")

# End ####
sessionInfo()