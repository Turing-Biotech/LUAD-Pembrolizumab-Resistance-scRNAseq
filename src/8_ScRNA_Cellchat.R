# Set ####
options(timeout = 7355608); getOption("timeout")
options(stringsAsFactors = FALSE)
getwd(); rm(list = ls()); gc()
set.seed(19981108)
if (!dir.exists("8_ScRNA_Cellchat")) {dir.create("8_ScRNA_Cellchat")}

# Library ####
library(Seurat); packageVersion("Seurat")
library(tidyverse); packageVersion("tidyverse")
library(patchwork); packageVersion("patchwork")
library(CellChat); packageVersion("CellChat")
library(ggsankey); packageVersion("ggsankey")
library(scales); packageVersion("scales")
library(grid); packageVersion("grid")

# Load ####
load(file = "1_ScRNA_SeuratData/Seurat_Celltype_Major.RData")
Seurat_Celltype_Major; head(Seurat_Celltype_Major@meta.data)
Idents(Seurat_Celltype_Major) <- Seurat_Celltype_Major@meta.data$Celltype_Major
levels(Idents(Seurat_Celltype_Major))
table(Seurat_Celltype_Major@meta.data$Group)
DefaultAssay(Seurat_Celltype_Major) <- "RNA"

# Seurat => CellChat ####
CellChat <- createCellChat(object = Seurat_Celltype_Major, group.by = "Celltype_Major", datatype = "RNA", assay = "RNA"); gc()
CellChat@DB <- CellChatDB.human

# Seurat => CellChat ####
Seurat_Nonresponse <- subset(Seurat_Celltype_Major, subset = Group == "Nonresponse")
Seurat_Response <- subset(Seurat_Celltype_Major, subset = Group == "Response")
CellChat_Nonresponse <- createCellChat(object = Seurat_Nonresponse, group.by = "Celltype_Major", datatype = "RNA", assay = "RNA")
CellChat_Response <- createCellChat(object = Seurat_Response, group.by = "Celltype_Major", datatype = "RNA", assay = "RNA")
levels(CellChat_Nonresponse@idents); levels(CellChat_Response@idents)
rm(Seurat_Celltype_Major, Seurat_Response, Seurat_Nonresponse); gc()
CellChat_Nonresponse@DB <- CellChatDB.human
CellChat_Response@DB  <- CellChatDB.human

# Run_CellChat_Pipeline ####
Run_CellChat_Pipeline <- function(CellChat_Object) {
  CellChat_Object <- subsetData(CellChat_Object)
  CellChat_Object <- identifyOverExpressedGenes(CellChat_Object)
  CellChat_Object <- identifyOverExpressedInteractions(CellChat_Object)
  CellChat_Object <- computeCommunProb(CellChat_Object)
  CellChat_Object <- filterCommunication(CellChat_Object, min.cells = 10)
  CellChat_Object <- computeCommunProbPathway(CellChat_Object)
  CellChat_Object <- aggregateNet(CellChat_Object)
  return(CellChat_Object)}
CellChat <- Run_CellChat_Pipeline(CellChat_Object = CellChat)
CellChat_Nonresponse <- Run_CellChat_Pipeline(CellChat_Object = CellChat_Nonresponse)
CellChat_Response <- Run_CellChat_Pipeline(CellChat_Object = CellChat_Response)

save(CellChat, file = "8_ScRNA_Cellchat/CellChat.RData")
save(CellChat_Nonresponse, file = "8_ScRNA_Cellchat/CellChat_Nonresponse.RData")
save(CellChat_Response, file = "8_ScRNA_Cellchat/CellChat_NO.RData")
load(file = "8_ScRNA_Cellchat/CellChat.RData")
load(file = "8_ScRNA_Cellchat/CellChat_Nonresponse.RData")
load(file = "8_ScRNA_Cellchat/CellChat_NO.RData")

pdf("8_ScRNA_Cellchat/Plot_Network_Count_All.pdf", width = 8, height = 10)
netVisual_circle(CellChat@net$count, 
                 weight.scale = TRUE, 
                 label.edge= FALSE, 
                 title.name = "Number of interactions"); dev.off()
pdf("8_ScRNA_Cellchat/Plot_Network_Weight_All.pdf", width = 8, height = 10)
netVisual_circle(CellChat@net$weight, 
                 weight.scale = TRUE, 
                 label.edge= FALSE, 
                 title.name = "Interaction weights/strength"); dev.off()
levels(CellChat@idents)

pdf("8_ScRNA_Cellchat/Plot_Network_Sender_Count_All.pdf", width = 8, height = 10)
netVisual_circle(CellChat@net$count, 
                 sources.use = c("T_Cell"), 
                 targets.use = NULL,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Number of interactions (Sender)"); dev.off() 
pdf("8_ScRNA_Cellchat/Plot_Network_Receiver_Count_All.pdf", width = 8, height = 10)
netVisual_circle(CellChat@net$count, 
                 sources.use = NULL, 
                 targets.use = c("T_Cell"),
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Number of interactions (Receiver)"); dev.off() 
pdf("8_ScRNA_Cellchat/Plot_Network_Sender_Weight_All.pdf", width = 8, height = 10)
netVisual_circle(CellChat@net$weight, 
                 sources.use = c("T_Cell"), 
                 targets.use = NULL,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Interaction weights/strength (Sender)"); dev.off()
pdf("8_ScRNA_Cellchat/Plot_Network_Receiver_Weight_All.pdf", width = 8, height = 10)
netVisual_circle(CellChat@net$weight, 
                 sources.use = NULL, 
                 targets.use = c("T_Cell"),
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Interaction weights/strength (Receiver)"); dev.off() 

# Plot_Bubble
levels(CellChat@idents)
Plot_Bubble_Source <- netVisual_bubble(object = CellChat,
                                       sources.use = c("T_Cell"),
                                       targets.use = c("B_Cell", "Monocyte", "NK_Cell", 
                                                       "Plasma_Cell"),
                                       angle.x = 45,
                                       remove.isolate = TRUE) + 
  guides(fill = guide_legend(override.aes = list(size = 10))) +
  labs(title = "Differential Interaction", y = "", x = "") + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.position = "right", 
        legend.key.size = unit(0.5, "in")); Plot_Bubble_Source
ggsave(plot = Plot_Bubble_Source, filename = "8_ScRNA_Cellchat/Plot_Bubble_Source_All.pdf", 
       width = 6, height = 10, dpi = 300, units = "in")

Plot_Bubble_Target <- netVisual_bubble(object = CellChat,
                                       sources.use = c("B_Cell", "Monocyte", "NK_Cell", 
                                                       "Plasma_Cell"),
                                       targets.use = c("T_Cell"),
                                       angle.x = 45,
                                       remove.isolate = TRUE) + 
  guides(fill = guide_legend(override.aes = list(size = 10))) +
  labs(title = "Differential Interaction", y = "", x = "") + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.position = "right", 
        legend.key.size = unit(0.5, "in")); Plot_Bubble_Target
ggsave(plot = Plot_Bubble_Target, filename = "8_ScRNA_Cellchat/Plot_Bubble_Target_All.pdf", 
       width = 6, height = 10, dpi = 300, units = "in")

slotNames(CellChat)
CellChat_LR <- subsetCommunication(
  object = CellChat)
head(CellChat_LR, n = 10)
Data_Sankey <- CellChat_LR %>%
  dplyr::filter(
    pval < 0.05, 
    source == "T_Cell", 
    target != "T_Cell") %>%  
  dplyr::select(source, target, interaction_name_2, ligand, receptor, pathway_name, annotation) %>%
  dplyr::distinct() %>%  
  ggsankey::make_long(source, target, interaction_name_2, ligand, receptor, pathway_name, annotation) %>%
  dplyr::mutate(node = factor(node))
Plot_Sankey <- ggplot(
  data = Data_Sankey, 
  mapping = aes(
    x = x, next_x = next_x, 
    node = node, next_node = next_node, 
    fill = node, label = node)) + 
  labs(title = "Receptor Ligands of Source (T Cell)") +
  geom_sankey(flow.alpha = 0.4, width = 0.12, smooth = 10) +
  geom_sankey_text(size = 4, type = "sankey", hjust = 0.5) +
  theme_void() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold")); Plot_Sankey
ggsave(plot = Plot_Sankey, filename = "8_ScRNA_Cellchat/Plot_Sankey_Source_All.pdf", 
       width = 15, height = 10, dpi = 300, units = "in")

# Merge CellChat ####
CellChat_List <- list(Nonresponse = CellChat_Nonresponse, Response = CellChat_Response)
CellChat_Group <- mergeCellChat(CellChat_List, add.names = names(CellChat_List))
print(names(CellChat_List)); CellChat_Group

# Plot_Bar_Count
Plot_Bar_Count <- compareInteractions(CellChat_Group, show.legend = TRUE, group = c("Nonresponse", "Response"), measure = "count") + 
  labs(title = "Barplot of Count Interaction", y = "Number of Inferres Interactions", x = "Group") + 
  scale_fill_manual(values = c("#FF6354FF", "#33CC33FF")) + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.key.size = unit(0.5, "in")); Plot_Bar_Count
ggsave(plot = Plot_Bar_Count, filename = "8_ScRNA_Cellchat/Plot_Bar_Count.pdf", 
       width = 6, height = 5, dpi = 300, units = "in")

# Plot_Bar_Weight
Plot_Bar_Weight <- compareInteractions(CellChat_Group, show.legend = TRUE, group = c("Nonresponse", "Response"), measure = "weight") +
  labs(title = "Barplot of Weight Interaction", y = "Interaction Strength", x = "Group") + 
  scale_fill_manual(values = c("#FF6354FF", "#33CC33FF")) + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.key.size = unit(0.5, "in")); Plot_Bar_Weight
ggsave(plot = Plot_Bar_Weight, filename = "8_ScRNA_Cellchat/Plot_Bar_Weight.pdf", 
       width = 6, height = 5, dpi = 300, units = "in")

# Plot_Bar_Relative
Plot_Bar_Relative <- rankNet(CellChat_Group, mode = "comparison", stacked = TRUE, do.stat = TRUE) +
  labs(title = "Relative Information Flow", y = "Signaling Percentage", x = "Signaling Pathway") + 
  scale_fill_manual(values = c("#33CC33FF", "#FF6354FF")) + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.position = "top", 
        legend.key.size = unit(0.5, "in")); Plot_Bar_Relative
ggsave(plot = Plot_Bar_Relative, filename = "8_ScRNA_Cellchat/Plot_Bar_Relative.pdf", 
       width = 6, height = 8, dpi = 300, units = "in")

# Plot_Network
print(names(CellChat_List)); CellChat_Group
CellChat_Group@net[["Nonresponse"]][["count"]]
CellChat_Group@net[["Response"]][["count"]]

pdf("8_ScRNA_Cellchat/Plot_Network_Count_Group.pdf", width = 8, height = 10)
netVisual_diffInteraction(object = CellChat_Group, 
                          color.edge = c("#FF6354FF", "#33CC33FF"), 
                          weight.scale = TRUE, 
                          measure = "count", 
                          shape = "circle", 
                          title.name = "Number of interactions", 
                          comparison = c("Response", "Nonresponse")); dev.off() 
pdf("8_ScRNA_Cellchat/Plot_Network_Weight_Group.pdf", width = 8, height = 10)
netVisual_diffInteraction(object = CellChat_Group, 
                          color.edge = c("#FF6354FF", "#33CC33FF"), 
                          weight.scale = TRUE, 
                          measure = "weight", 
                          shape = "circle", 
                          title.name = "Interaction weights/strength", 
                          comparison = c("Response", "Nonresponse")); dev.off() 
pdf("8_ScRNA_Cellchat/Plot_Network_Sender_Count_Group.pdf", width = 8, height = 10)
netVisual_diffInteraction(object = CellChat_Group, 
                          color.edge = c("#FF6354FF", "#33CC33FF"), 
                          weight.scale = TRUE, 
                          measure = "count", 
                          shape = "circle", 
                          sources.use = c("T_Cell"), 
                          targets.use = c("B_Cell", "Monocyte", "NK_Cell", 
                                          "Plasma_Cell", "Platelet"), 
                          title.name = "Number of interactions (Sender)", 
                          comparison = c("Response", "Nonresponse")); dev.off() 
pdf("8_ScRNA_Cellchat/Plot_Network_Sender_Weight_Group.pdf", width = 8, height = 10)
netVisual_diffInteraction(object = CellChat_Group, 
                          color.edge = c("#FF6354FF", "#33CC33FF"), 
                          weight.scale = TRUE, 
                          measure = "weight", 
                          shape = "circle", 
                          sources.use = c("T_Cell"), 
                          targets.use = c("B_Cell", "Monocyte", "NK_Cell", 
                                          "Plasma_Cell", "Platelet"), 
                          title.name = "Interaction weights/strength (Sender)", 
                          comparison = c("Response", "Nonresponse")); dev.off() 
pdf("8_ScRNA_Cellchat/Plot_Network_Receiver_Count_Group.pdf", width = 8, height = 10)
netVisual_diffInteraction(object = CellChat_Group, 
                          color.edge = c("#FF6354FF", "#33CC33FF"), 
                          weight.scale = TRUE, 
                          measure = "count", 
                          shape = "circle", 
                          sources.use = c("B_Cell", "Monocyte", "NK_Cell", 
                                          "Plasma_Cell", "Platelet"), 
                          targets.use = c("T_Cell"), 
                          title.name = "Number of interactions (Receiver)", 
                          comparison = c("Response", "Nonresponse")); dev.off() 
pdf("8_ScRNA_Cellchat/Plot_Network_Receiver_Weight_Group.pdf", width = 8, height = 10)
netVisual_diffInteraction(object = CellChat_Group, 
                          color.edge = c("#FF6354FF", "#33CC33FF"), 
                          weight.scale = TRUE, 
                          measure = "weight", 
                          shape = "circle", 
                          sources.use = c("B_Cell", "Monocyte", "NK_Cell", 
                                          "Plasma_Cell", "Platelet"), 
                          targets.use = c("T_Cell"), 
                          title.name = "Interaction weights/strength (Receiver)", 
                          comparison = c("Response", "Nonresponse")); dev.off() 

# Plot_Bubble
Plot_Bubble_Source_Group <- netVisual_bubble(object = CellChat_Group,
                                           sources.use = c("T_Cell"),
                                           targets.use = c("B_Cell", "Monocyte", "NK_Cell", 
                                                           "Plasma_Cell"),
                                           comparison = c(1, 2),
                                           angle.x = 45,
                                           remove.isolate = TRUE) + 
  guides(fill = guide_legend(override.aes = list(size = 10))) +
  labs(title = "Differential Interaction (Nonresponse vs Response)", y = "", x = "") + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.position = "right", 
        legend.key.size = unit(0.5, "in")); Plot_Bubble_Source_Group
ggsave(plot = Plot_Bubble_Source_Group, filename = "8_ScRNA_Cellchat/Plot_Bubble_Source_Group.pdf", 
       width = 6, height = 10, dpi = 300, units = "in")

Plot_Bubble_Target_Group <- netVisual_bubble(object = CellChat_Group,
                                           sources.use = c("B_Cell", "Endothelial_Cell", "Endothelial_Cell_Lymphatic", 
                                                           "Epithelial_Cell", "Myeloid_Cell", "Smooth_Muscle_Cell", 
                                                           "T_Cell", "Granulosa_Cell", "Stromal_Cell"),
                                           targets.use = c("T_Cell"),
                                           comparison = c(1, 2),
                                           angle.x = 45,
                                           remove.isolate = TRUE) + 
  guides(fill = guide_legend(override.aes = list(size = 10))) +
  labs(title = "Differential Interaction (Nonresponse vs Response)", y = "", x = "") + 
  theme(axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.position = "right", 
        legend.key.size = unit(0.5, "in")); Plot_Bubble_Target_Group
ggsave(plot = Plot_Bubble_Target_Group, filename = "8_ScRNA_Cellchat/Plot_Bubble_Target_Group.pdf", 
       width = 5, height = 10, dpi = 300, units = "in")

# File Tree ####
fs::dir_tree("8_ScRNA_Cellchat/", max_depth = 2)

# End ####
sessionInfo()
