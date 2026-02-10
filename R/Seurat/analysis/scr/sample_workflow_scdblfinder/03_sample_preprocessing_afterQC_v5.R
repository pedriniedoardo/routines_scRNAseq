# AIM ---------------------------------------------------------------------
# Apply the filter defined in the LUT_QC_after.csv file
# This is either a testing of the filtering results, or the actual filtering choice to go forward
# this is the processing with scDblFinder

# renv integration --------------------------------------------------------

# to load the packages
source(".Rprofile")

# in the config specify the following
# renv_library_path: "renv/library/linux-rocky-9.5/R-4.5/x86_64-conda-linux-gnu"

# libraries ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(scDblFinder)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1000 * 1024^2)

# define the filtering parameters -----------------------------------------
# read in the file with the filtering parameters
LUT_QC_after <- read_csv("data/LUT_QC_after.csv")
label <- "afterQC"

# read in the data --------------------------------------------------------
# location of all the raw matrices
list_sobj <- readRDS("out/object/analysis_R45/02_list_sobj.rds")

# check the dimensions before filtering
lapply(list_sobj,function(x){dim(x)})

# processing --------------------------------------------------------------

# scobj <- list_sobj[[1]]
# nm <- "s1"
list_sobj2 <- pmap(list(list_sobj,names(list_sobj)),function(scobj,nm){
  # track the progress
  print(nm)
  
  # define the filteirng parameters
  # featureLow_thr <- 500
  # featureHigh_thr <- 5000
  # mito_thr <- 10
  # label <- paste(featureLow_thr,featureHigh_thr,mito_thr,"V5",sep = "_")
  featureLow_thr <- LUT_QC_after %>% filter(sample_name == nm) %>% pull(featureLow_thr)
  featureHigh_thr <- LUT_QC_after %>% filter(sample_name == nm) %>% pull(featureHigh_thr)
  mito_thr <- LUT_QC_after %>% filter(sample_name == nm) %>% pull(mito_thr)
  
  # define the label based on the filtering parameters
  label <- paste(featureLow_thr,featureHigh_thr,mito_thr,"V5",sep = "_")
  
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
  # only remove low coverage (I have already used a filter of 200 features per barcode in the loading step) cells and proprocess to generata clusters needed for the doublet identification
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
  # set seed for reproducibility
  set.seed(123)
  sce_scobj <- scDblFinder(GetAssayData(scobj, layer="counts"), clusters=Idents(scobj))
  
  # port the resulting scores back to the Seurat object:
  scobj$scDblFinder.score <- sce_scobj$scDblFinder.score
  scobj$scDblFinder.class <- sce_scobj$scDblFinder.class
  
  # perform the filtering based on the fixed threshold defined
  # perform the filtering based on the doublets assignamente
  scobj_filter <- subset(scobj, subset = discard_threshold == 0 & scDblFinder.class == "singlet")
  
  # preprocess data after filtering 
  scobj_filter <- scobj_filter %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters() %>%
    RunUMAP(dims = 1:30) 
  
  return(scobj_filter)

}) %>%
  setNames(names(list_sobj))

# check the dimenstions after filtering
lapply(list_sobj2,function(x){dim(x)})

# save the list -----------------------------------------------------------
#
saveRDS(list_sobj2,paste0("out/object/analysis_R45/03_list_sobj_filtered_",label,".rds"))
