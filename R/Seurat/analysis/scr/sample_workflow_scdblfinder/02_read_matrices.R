# AIM ---------------------------------------------------------------------
# generate the sample rds for all the objects
# This is either a testing of the filtering results, or the actual filtering choice to go forward
# this is the processing with scDblFinder

# libraries ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(scuttle)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
# location of all the raw matrices
folder <- "data/raw/cellranger8/SoupX/"
id_sample <- dir(folder)

# load the LUT of the dataset
LUT <- read_csv("data/LUT_samples.csv")

# sample processing -------------------------------------------------------
# id <- "s1"
list_sobj <- lapply(id_sample,function(id){
  # keep the tracking
  print(id)
  
  # identity the input files
  file_id <- file.path(folder,id)
  
  # read in the matrix
  data <- Read10X(data.dir = file_id)
  
  # crete the object
  datasc <- CreateSeuratObject(counts = data, project = id,
                               # potentially parametrize the values below
                               # remove low expressed features
                               # in this case keep all the features before the integration
                               min.cells = 0,
                               # remove low content barcodes
                               min.features = 200)
  
  # add the QC metrics to the metadata
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
  datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
  
  # add more covariates from the LUT
  datasc$pathology <- LUT %>%
    filter(sample == id) %>%
    pull(pathology)
  
  datasc$location <- LUT %>%
    filter(sample == id) %>%
    pull(location)
  
  datasc$tissue <- LUT %>%
    filter(sample == id) %>%
    pull(tissue)
  
  # -------------------------------------------------------------------------
  # add the filtering variable based on the adaptive threshold multivalue
  stats <- cbind(log10(datasc@meta.data$nCount_RNA),
                 log10(datasc@meta.data$nFeature_RNA),
                 datasc@meta.data$percent.mt)
  
  outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
  multi.outlier <- isOutlier(outlying, type = "higher")
  datasc$discard_multi <- as.vector(multi.outlier)
  
  # add the filtering variable based on the adaptive threshold single values
  high_QC_mito <- isOutlier(datasc@meta.data$percent.mt, type="high", log=TRUE)
  QC_features <- isOutlier(datasc@meta.data$nFeature_RNA, type="both", log=TRUE)
  
  datasc$discard_single <- high_QC_mito | QC_features
  
  return(datasc)
  # saveRDS(datasc,file = paste0("../../out/object/02_",id,",rds"))
  # write_tsv(datasc@meta.data %>% rownames_to_column("barcodes"),paste0("../../out/table/02_",id,".tsv"))
}) %>%
  setNames(id_sample)

# save the list -----------------------------------------------------------
#
saveRDS(list_sobj,"out/object/analysis_R45/02_list_sobj.rds")
