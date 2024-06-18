rm(list = ls()) # remove environment variables
library(tidyverse)
library(Seurat)
library(harmony)

# import data
sce.mer <- readRDS("/data/whu/home/whr_brain_project_final/results/base/hippocampus_harmony_final.rds")
table(sce.mer$celltype)
table(sce.mer$orig.ident)

Resolution = 0.5
Dir <- "/data/whu/home/whr_brain_project_final/"
dir.create(paste0(Dir, "figures/figure_3/subcluster_", Resolution), recursive = T)
dir.create(paste0(Dir, "results/figure_3/subcluster_", Resolution), recursive = T)

for(Cell in  c("ExN CA", "ExN DG", "InN", "MG", "ASC", "ODC", "OPC")){
  print(paste0(Cell, ": ", Resolution))
  
  figure_folder <- paste0(Dir, "figures/figure_3/subcluster_", Resolution, "/", Cell)
  if(!dir.exists(figure_folder)){
    dir.create(figure_folder)
  }
  
  result_folder <- paste0(Dir, "results/figure_3/subcluster_", Resolution, "/", Cell)
  if(!dir.exists(result_folder)){
    dir.create(result_folder)
  }
  
  ## extract data of cell types
  sce.sub <- subset(sce.mer, celltype == Cell)
  if(Cell %in% c("MG", "ASC")){
    sce.sub <- subset(sce.sub, orig.ident != "30H")} # remove samples without > 10 cells
  df <- CreateSeuratObject(counts = sce.sub@assays$RNA@counts,  meta.data = sce.sub@meta.data, min.cells = 3, min.features = 200)
  
  ## RunHarmony and clustering 
  df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
  df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)
  df <- ScaleData(df, features = rownames(df))
  df <- RunPCA(df, features = VariableFeatures(object = df))
  df <- df %>%
    RunHarmony("orig.ident", plot_convergence = FALSE)
  df <- df %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = Resolution) %>%
    RunUMAP(reduction = "harmony", dims = 1:30) %>%
    identity()
  
  # export seurat object
  setwd(result_folder)
  saveRDS(df, file = paste0(Cell, ".rds"))
}
