# Set ####
options(timeout = 7355608); getOption("timeout")
options(stringsAsFactors = FALSE)
getwd(); rm(list = ls()); gc()
if (!dir.exists("1_ScRNA_SeuratData")) {dir.create("1_ScRNA_SeuratData")}

# Library ####
library(Seurat); packageVersion("Seurat")
library(tidyverse); packageVersion("tidyverse")
library(bruceR); packageVersion("bruceR")

# Seurat ####
File <- list.files(path = "0_ScRNA_RawData", full.names = TRUE)
Seurat_List <- list()
for (Sample in File) {
  Data <- Read10X(data.dir = Sample, 
                  gene.column = 2, # features.tsv.gz
                  cell.column = 1, # barcodes.tsv.gz
                  unique.features = TRUE) 
  Name <- basename(Sample)
  Seurat <- CreateSeuratObject(counts = Data, 
                               project = Name,
                               min.cells = 3, 
                               min.features = 200)
  Seurat_List[[Name]] <- Seurat
  rm(Data, Seurat); gc()}
Seurat_List; names(Seurat_List)

# Merge ####
Seurat <- merge(x = Seurat_List[[1]], 
                y = Seurat_List[-1], 
                add.cell.ids = names(Seurat_List), 
                project = "LUAD") # [!]
Seurat; dim(Seurat) # 29453 Feature; 106385 Cell
head(Seurat@meta.data, n = 10); table(Seurat@meta.data$orig.ident)

# Summary ####
Seurat_Summary <- Seurat@meta.data %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(Cell = n()) %>% 
  dplyr::rename(Sample = orig.ident) %>% 
  dplyr::mutate(Gender = dplyr::case_when(
    Sample %in% "BH1" ~ "Nonresponse_1",
    Sample %in% "BH2" ~ "Nonresponse_2",
    Sample %in% "BH4" ~ "Nonresponse_3",
    Sample %in% "BH5" ~ "Response_1",
    Sample %in% "BH7" ~ "Response_2",
    Sample %in% "BH16" ~ "Response_3",
    Sample %in% "FH1" ~ "Nonresponse_4",
    TRUE ~ "Unknown")); Seurat_Summary
save(Seurat_Summary, file = "1_ScRNA_SeuratData/Seurat_Summary.RData")
print_table(Seurat_Summary)

# Meta_Data Tidy ####
# Sample
table(Seurat@meta.data$orig.ident)
Seurat@meta.data <- Seurat@meta.data %>% 
  dplyr::mutate(Sample = dplyr::case_when(
    Seurat@meta.data$orig.ident %in% "BH1" ~ "Nonresponse_1",
    Seurat@meta.data$orig.ident %in% "BH2" ~ "Nonresponse_2",
    Seurat@meta.data$orig.ident %in% "BH4" ~ "Nonresponse_3",
    Seurat@meta.data$orig.ident %in% "BH5" ~ "Response_1",
    Seurat@meta.data$orig.ident %in% "BH7" ~ "Response_2",
    Seurat@meta.data$orig.ident %in% "BH16" ~ "Response_3",
    Seurat@meta.data$orig.ident %in% "FH1" ~ "Nonresponse_4",
    TRUE ~ "Unknown"))
head(Seurat@meta.data, n = 10); table(Seurat@meta.data$Sample)

# Group
Seurat@meta.data <- Seurat@meta.data %>% 
  dplyr::mutate(Group = dplyr::case_when(
    Seurat@meta.data$Sample %in% c("Nonresponse_1", "Nonresponse_2", "Nonresponse_3", "Nonresponse_4") ~ "Nonresponse",
    Seurat@meta.data$Sample %in% c("Response_1", "Response_2", "Response_3") ~ "Response",
    TRUE ~ "Unknown"))
head(Seurat@meta.data, n = 10); table(Seurat@meta.data$Group)

# Save
save(Seurat_List, file = "1_ScRNA_SeuratData/Seurat_List.RData")
save(Seurat, file = "1_ScRNA_SeuratData/Seurat.RData")

# End ####
sessionInfo()
