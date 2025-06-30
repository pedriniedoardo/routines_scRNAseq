# AIM ---------------------------------------------------------------------
# for this step the version of Seurat is not important.
# transfer the minimal number of files from the source directory to the project folder. Make sure the folder analysis/ filtered_feature_bc_matrix/ and raw_feature_bc_matrix/ are present in the project folder.
# run SoupX and generate the cleaned raw matrices.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(SoupX)
library(DropletUtils)

# run the soupx processing ------------------------------------------------
id_sample <- dir("../data/test_introns/cellranger7_out/")

lapply(id_sample,function(x){
  # track the progress
  print(x)
  # define the location of the output of cellranger
  file <- paste0("../data/test_introns/cellranger7_out/",x)
  # Load data and estimate soup profile
  sc <- load10X(file)
  # Estimate rho
  sc <- autoEstCont(sc)
  # Clean the data
  out <- adjustCounts(sc,roundToInt = T)
  # save the data
  DropletUtils:::write10xCounts(paste0("../data/test_introns/SoupX_default_cellranger7/",x), out)
})
