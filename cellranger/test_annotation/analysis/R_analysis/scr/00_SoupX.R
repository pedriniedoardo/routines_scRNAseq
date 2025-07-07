# AIM ---------------------------------------------------------------------
# for this step the version of Seurat is not important.
# transfer the minimal number of files from the source directory to the project folder. Make sure the folder analysis/ filtered_feature_bc_matrix/ and raw_feature_bc_matrix/ are present in the project folder.
# run SoupX and generate the cleaned raw matrices.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(SoupX)
library(DropletUtils)

# run the soupx processing ------------------------------------------------
# define the location of the output of cellranger
cellranger_out <- "../../out/cellranger901/"
SopuX_out <- "../../data/SopuX_out/"

id_sample <- dir(cellranger_out) %>%
  # str_subset(pattern = "CT|DS") %>%
  str_subset(pattern = ".log",negate = T)

# x <- "SacsKO5mo_MG"
lapply(id_sample,function(x){
  # track the progress
  print(x)
  # define the location of the output of cellranger
  file <- paste0(cellranger_out,x,"/outs/")
  # Load data and estimate soup profile
  sc <- load10X(file)
  # Estimate rho
  sc <- autoEstCont(sc)
  # Clean the data
  out <- adjustCounts(sc,roundToInt = T)
  # save the data
  DropletUtils:::write10xCounts(paste0(SopuX_out,x), out)
})
