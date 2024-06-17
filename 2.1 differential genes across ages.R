rm(list = ls()) # remove environment variables
library(tidyverse)
library(monocle)

cell_types <- c("InN", "ExN DG", "ASC", "ODC", "OPC", "EC", "PC", "VLMC", "EP")
core <- 1

for(Cell in cell_types){
  print(Cell)
  
  if(Cell != "VLMC"){
    cluster <- 25
  }else{
    cluster <- 5
  }
  
  Dir <- "/data/whu/home/whr_brain_project_final/"
  figure_folder <- paste0(Dir, "figures/figure_2/age_deg_adjusted_both/", Cell)
  if(!dir.exists(figure_folder)){
    dir.create(figure_folder, recursive = T)
  }
  
  result_folder <- paste0(Dir, "results/figure_2/age_deg_adjusted_both/", Cell)
  if(!dir.exists(result_folder)){
    dir.create(result_folder, recursive = T)
  }
  
  df <- readRDS(paste0("/data/whu/home/whr_brain_project_final/results/figure_3/subcluster_0.5/", Cell, "/", Cell, ".rds"))
  
  # add PMI group
  df$PMI_group <- NA
  df$PMI_group[df$PMI <= 9] <- "4-9 hour"
  df$PMI_group[df$PMI > 9 & df$PMI <= 15] <- "9-15 hour"
  df$PMI_group[df$PMI > 15] <- "20 hour"
  df$PMI_group <- factor(df$PMI_group, c("4-9 hour", "9-15 hour", "20 hour"))
  
  exp <- df@assays$RNA@data
  
  meta <- df@meta.data
  pd <- new('AnnotatedDataFrame', data = meta)
  fData <- data.frame(gene_short_name = row.names(exp), row.names = row.names(exp))
  fd <- new('AnnotatedDataFrame', data = fData)
  monocle_cds <- newCellDataSet(exp,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = uninormal())
  
  monocle_cds$Pseudotime <- monocle_cds$age
  time_diff <- differentialGeneTest(monocle_cds, cores = core, 
                                    fullModelFormulaStr = "~sm.ns(Pseudotime) + gender + PMI_group",
                                    reducedModelFormulaStr = "~gender + PMI_group")
  time_diff$pval_FDR <- p.adjust(time_diff$pval, method = "fdr")
  time_diff$qval_FDR <- p.adjust(time_diff$qval, method = "fdr")
  
  setwd(result_folder)
  write.table(time_diff, paste0(Cell, "_age_deg.txt"), row.names = F, sep="\t", quote = F)
  saveRDS(monocle_cds, file = paste0(Cell, "_monocle.rds"))
  
  mat_col <- data.frame(Pseudotime = seq(min(monocle_cds$Pseudotime), max(monocle_cds$Pseudotime), length.out = 100))
  colnames(mat_col) <- "Age"
  newdata <- mat_col
  colnames(newdata) <- "Pseudotime"
  
  time_diff_sig <- subset(time_diff, pval_FDR < 0.01 & qval_FDR < 0.01) # set adjusted p value and q value both < 0.01
  p1 <- plot_pseudotime_heatmap(monocle_cds[time_diff_sig$gene_short_name,],
                                num_clusters = cluster,
                                cluster_rows = T,
                                add_annotation_col = mat_col,
                                show_rownames = F, 
                                return_heatmap = T,
                                hmcols = colorRampPalette(colors = c("blue", "white", "red"))(100),
                                cores = core)
  setwd(figure_folder)
  ggsave(paste0(Cell, "_age_sig_gene_heatmap.pdf"), p1)
  
  
  genSmoothCurves_sig <- genSmoothCurves(monocle_cds[time_diff_sig$gene_short_name,],
                                         new_data = newdata,
                                         cores = core)
  
  clusters <- data.frame(cutree(p1$tree_row, k = cluster))
  colnames(clusters) <- "gene_clusters"
  clusters <- data.frame(gene = rownames(clusters), clusters)
  
  setwd(result_folder)
  write.table(clusters, paste0(Cell, "_age_sig_gene_cluster.txt"), row.names = F, sep="\t", quote = F)
}

