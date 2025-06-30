# AIM ---------------------------------------------------------------------
# Run the preprocessing of an individula sample.
# The aim is to save the pre-QC filtering object and the pre-QC metadata table.

# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(scuttle)
# library(homologene)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# locate the input --------------------------------------------------------
# define the project folder
project_folder <- "../data/"

# define the sample id
in_id_sample <- "connect_5k_pbmc_NGSC3_ch1_gex_1"

# define the data to be loaded
in_id_data <- "outs/filtered_feature_bc_matrix"

# define ogganisim
in_id_org <- "9606"

# potentially load the LUT for the sample
# LUT <- read_csv("../data/test_introns/LUT_samples.csv")

# run the processing ------------------------------------------------------
# read in the matrix
data <- Read10X(data.dir = file.path(project_folder,in_id_sample,in_id_data))
  
# crete the object
datasc <- CreateSeuratObject(counts = data, project = in_id_sample,
                             # potentially parametrize the values below
                             min.cells = 20,
                             min.features = 200)

# add the mefadata accordin to the organism
if(in_id_org == "9606"){
  # add the metadata
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
  datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
} else if(in_id_org == "10090"){
  # add the metadata
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^Mt-")
  datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")
  datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^Hb[^(p)]")
} else {
  # decide whether to run something or nothing
}


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
# -------------------------------------------------------------------------

# saving ------------------------------------------------------------------
# save the object pre application of the QC filters
out_id_object <- paste0("01_",in_id_sample,"_obj_preQC.rds")
out_id_meta <- paste0("01_",in_id_sample,"_meta_preQC.tsv")

saveRDS(datasc,file.path("../out/object",out_id_object))
write_tsv(datasc@meta.data %>% rownames_to_column("barcodes"),file.path("../out/table",out_id_meta))
