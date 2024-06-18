rm(list = ls()) # remove environment variables
library(tidyverse)
library(Seurat)
library(harmony)

######## Create Seurat objects ########

# get sample list
dir_matrix <- "/data/whu/home/whr_brain_project_final/results/base/splited_strainedCounts"
samples <- list.files(dir_matrix)

samList <- lapply(samples, function(sp){
  CreateSeuratObject(counts = Read10X_h5(file.path(dir_matrix, sp)),
                     project = gsub("\\..*", "", sp),
                     min.cells = 3,
                     min.features = 200)
})
names(samList) <- gsub("\\..*", "", samples)

######## Filter cells ########

dir_scrublet="/data/whu/home/whr_brain_project_final/data/adjusted_scrublet/"

rt_summary <- data.frame()
for(i in names(samList)){
  samList[[i]] <- PercentageFeatureSet(object = samList[[i]],  pattern = "^MT-", col.name = "percent.mt")
  cell_before <- length(colnames(samList[[i]]))
  
  scrublet_result <- read.table(paste0(dir_scrublet, i,"_scrublet_results.tsv"), header = T, sep = "\t", check.names = F)
  scrublet_result <- scrublet_result[match(colnames(samList[[i]]), scrublet_result$Barcode),]
  
  if(identical(scrublet_result$Barcode, colnames(samList[[i]]))){
    samList[[i]]$doublet <- scrublet_result$scrublet_DropletType
    samList[[i]]$doublet_score <- scrublet_result$scrublet_Scores
  }else{
    print("Wrong sequence")
  }
  
  doublet_number <- sum(samList[[i]]$doublet == "doublet")
  
  samList[[i]] <- subset(samList[[i]], nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & doublet == "singlet")
  
  cell_after <- length(colnames(samList[[i]]))
  
  # summary
  rt_summary_i <- data.frame(Sample = i, Before_Filter = cell_before, After_Filter = cell_after, Doublet_number = doublet_number)
  rt_summary <- rbind(rt_summary, rt_summary_i)
  
  cat("Sample: ", i, "; Before filter:", cell_before, "cells; ", "After filter:", cell_after,"cells; Doublet number: ", doublet_number, "\n")
  
}

rt_summary$Doublet_rate <- rt_summary$Doublet_number*100/rt_summary$Before_Filter

# merge data
df <- samList[[1]]
for (i in 2:length(samList)) {
  print(i)
  df <- merge(df, y = samList[[i]])
}

## RunHarmony and clustering 
df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)
df <- ScaleData(df, features = rownames(df))
df <- RunPCA(df, features = VariableFeatures(object = df))
df <- df %>%
  RunHarmony("orig.ident", plot_convergence = FALSE)
df <- df %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  identity()

# export seurat object
setwd("/data/whu/home/whr_brain_project_final/results/base")
saveRDS(df, file = "hippocampus_harmony.rds")
write.table(rt_summary, file = "summary_running.txt", sep = "\t", row.names = F, quote = F)
