# AIM ---------------------------------------------------------------------
# This script is run to explore the dataset individually and define a consensus threshold for the mito and reads counts
# This is run on the post SoupX data correction.

# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# READ IN DATA ------------------------------------------------------------
folder <- "../../out/cellranger901/"
id_sample <- dir(folder) %>%
  str_subset(pattern = "test_neuron")

# data.frame(sample = id_sample) %>%
#   separate(col = "sample",into = c("treat","sample_id"),sep = "_",remove = F) %>%
#   write_csv("../data/LUT_samples.csv")

# load the LUT
LUT <- read_csv("../../data/LUT_sample.csv")

# do the preprocessing over all the dataset and save the objects
# x <- "test_neuron_auto"
list_datasc <- lapply(id_sample,function(x){
  # to track the processing of the progress of the lapply
  print(x)
  
  # read in the matrix
  data <- Read10X(data.dir = paste0(folder,x,"/outs/filtered_feature_bc_matrix/"))
  
  # crete the object
  datasc <- CreateSeuratObject(counts = data, project = LUT %>%
                                 filter(sample_id == x) %>%
                                 pull(sample_name), min.cells = 20, min.features = 200)
  
  # add the metadata
  # datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
  # datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^mt-")
  datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")
  # add also the percentage of globin. in this dataset it is not meaningful as there is no blood
  # datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
  datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^Hb[^(p)]")
  
  # label the cells based on the mt reads content
  datasc$mt_bin <- datasc@meta.data %>%
    mutate(test = case_when(percent.mt < 1~"low",
                            percent.mt < 10~"mid",
                            T ~ "high")) %>%
    pull(test)
  
  # datasc$treat <- LUT %>%
  #   filter(sample_id == x) %>%
  #   pull(status)
  # 
  # datasc$ID <- LUT %>%
  #   filter(sample_id == x) %>%
  #   pull(sample_id)
  
  # add the filtering variable based on the fixed threshold. at this poit add some generic thresholds
  datasc$discard_threshold <- datasc@meta.data %>%
    # mutate(test = percent.mt > mito_thr | nFeature_RNA < featureLow_thr | nFeature_RNA > featureHigh_thr) %>% 
    mutate(test = percent.mt > 15 | nFeature_RNA < 1000 | nFeature_RNA > 6000) %>% 
    pull(test)
  
  # add the filtering variable based on the adaptive threshold, multivariable
  # library(robustbase)
  # library(scater)
  stats <- cbind(log10(datasc@meta.data$nCount_RNA),
                 log10(datasc@meta.data$nFeature_RNA),
                 datasc@meta.data$percent.mt)
  
  outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
  multi.outlier <- isOutlier(outlying, type = "higher")
  datasc$discard_multi <- as.vector(multi.outlier)
  
  # add the filtering variable based on the adaptive threshold single variable
  high_QC_mito <- isOutlier(datasc@meta.data$percent.mt, type="high", log=TRUE)
  QC_features <- isOutlier(datasc@meta.data$nFeature_RNA, type="both", log=TRUE)
  
  datasc$discard_single <- high_QC_mito | QC_features
  
  return(datasc)
}) %>%
  setNames(id_sample)

# confirm the class of the objects ----------------------------------------
lapply(list_datasc, function(x){
  class(x@assays$RNA)
})

lapply(list_datasc, function(x){
  dim(x)
})

all.equal(list_datasc$test_neuron_auto@assays$RNA$counts,list_datasc$test_neuron_token@assays$RNA$counts)
