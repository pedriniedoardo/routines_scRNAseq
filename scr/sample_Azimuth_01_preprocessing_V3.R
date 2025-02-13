# AIM ---------------------------------------------------------------------
# sample processing of Azimuth processing for reference labelling
# this part produce the file needed for the actula label transfer 

# libraries ---------------------------------------------------------------
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(tidyverse)
library(ComplexHeatmap)
library(cowplot)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v3")

# read in the dataset -----------------------------------------------------
# this is the query dataset
data.combined <- readRDS(file = "../out/test_introns/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15_V3.rds")

DimPlot(data.combined)

# generate the diet object and save it
# save the coordiantes of the UMAP
data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode") %>%
  write_tsv("../out/test_introns/table/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15_V3_coordUMAP.tsv")

DefaultAssay(data.combined) <- "RNA"
object_diet <- DietSeurat(object = data.combined, assays = "RNA")
saveRDS(object_diet,"../out/test_introns/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15_V3_diet.rds")

# available datasets ------------------------------------------------------
# pick a reference
available_data <- AvailableData()

# install the picked reference dataset
# InstallData("humancortexref")

# if not working it is possible to install manually from the tar.gz file
# download.file(url = "http://seurat.nygenome.org/src/contrib/humancortexref.SeuratData_1.0.0.tar.gz",
#               destfile = "../data/azimuth/humancortexref.SeuratData_1.0.0.tar.gz")
# install.packages("../data/azimuth/humancortexref.SeuratData_1.0.0.tar.gz")

# load the reference dataset
reference <- LoadReference(path = "renv/library/R-4.3/x86_64-pc-linux-gnu/humancortexref.SeuratData/azimuth/")

# save the reference umap
df_point <- reference$map@reductions$refUMAP@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")
write_tsv(df_point,"../out/test_introns/table/azimuth_humancortex_CoordUMAP.tsv")

# save the meta of the ref
df_meta <- reference$map@meta.data %>%
  data.frame() %>%
  rownames_to_column("barcode")
write_tsv(df_meta,"../out/test_introns/table/azimuth_humancortex_Metadata.tsv")
