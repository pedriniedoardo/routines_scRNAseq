# AIM ---------------------------------------------------------------------
# Apply the filter stated in the label of the file.
# This is either a testing of the filtering results, or the actual filtering choice to go forward

# libraries ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
# library(DoubletFinder)
library(scDblFinder)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1000 * 1024^2)

# define the filtering parameters -----------------------------------------
featureLow_thr <- 1000
featureHigh_thr <- 6000
mito_thr <- 15
# "01000_06000_15_V5"
# build the label
label <- paste(featureLow_thr,featureHigh_thr,mito_thr,"V5",sep = "_")

# read in the data --------------------------------------------------------
# define the sample id
in_id_sample <- "connect_5k_pbmc_NGSC3_ch1_gex_2"

# id of the object
id_scobj <- "../out/object/01_connect_5k_pbmc_NGSC3_ch1_gex_2_obj_preQC.rds"

# load the LUT of the dataset
# LUT <- read_csv("../data/test_introns/LUT_samples.csv")

# run the processing ------------------------------------------------------
scobj <- readRDS(id_scobj)

# add the filtering label
scobj$label <- label

# add the filtering variable based on the fixed threshold
scobj$discard_threshold <- scobj@meta.data %>%
  mutate(test = percent.mt > mito_thr | nFeature_RNA < featureLow_thr | nFeature_RNA > featureHigh_thr) %>% 
  pull(test)

# preprocess the dataset before the doublet identification as recommended in:
# https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
# 1.5.11
# according to the documentation the doublet identification should be run before the fine filtering
# remove low coverage cells and proprocess to generata clusters needed for the doublet identification
scobj <- subset(scobj,subset = nCount_RNA > 500) %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  # do I run the UMAP ? I do not need it for the doublet identification, but can be useful in case someone wants to explore an individual sample
  RunUMAP(dims = 1:30)

# run scDblFinder after filtering the low coverage cells
sce_scobj <- scDblFinder(GetAssayData(scobj, layer="counts"), clusters=Idents(scobj))

# port the resulting scores back to the Seurat object:
scobj$scDblFinder.score <- sce_scobj$scDblFinder.score
scobj$scDblFinder.class <- sce_scobj$scDblFinder.class

# perform the filtering based on the fixed threshold defined
scobj_filter <- subset(scobj, subset = discard_threshold == 0)

# -------------------------------------------------------------------------
# preprocess data after filtering ?
scobj_filter <- scobj_filter %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30) 
# -------------------------------------------------------------------------

# saving ------------------------------------------------------------------
# save the object after application of the QC filters
out_id_object <- paste0("02_",in_id_sample,"_obj_postQC.rds")
out_id_meta <- paste0("02_",in_id_sample,"_meta_postQC.tsv")

saveRDS(scobj_filter,file.path("../out/object",out_id_object))
write_tsv(scobj_filter@meta.data %>% rownames_to_column("barcodes"),file.path("../out/table",out_id_meta))
