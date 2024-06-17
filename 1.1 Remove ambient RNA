rm(list = ls()) # remove environment variables
library(tidyverse)
library(SoupX)

# build folders
result_folder <- "/data/whu/home/whr_brain_project_SoupX_2_splited_adjusted/results/base"
if(!dir.exists(result_folder)){
  dir.create(result_folder, recursive = T)
}

# get sample list
dir_matrix <- "/data/whu/home/whr_brain_project_2/data/matrix"
samples <- list.files(dir_matrix)

for(i in samples){
  print(i)
  
  # Load data and estimate soup profile
  sc = load10X(file.path(dir_matrix, i))
  
  # Estimate rho
  sc = autoEstCont(sc)
  
  # Clean the data
  out = adjustCounts(sc)
  
  # export
  saveRDS(sc, file.path("/data/whu/home/whr_brain_project_SoupX_2_splited_adjusted/results/base/SoupX", paste0(i, ".rds")))
  DropletUtils:::write10xCounts(file.path("/data/whu/home/whr_brain_project_SoupX_2_splited_adjusted/results/base/splited_strainedCounts", paste0(i, ".h5")), out, version="3", type="HDF5")
}
