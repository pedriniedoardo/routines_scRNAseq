# AIM ---------------------------------------------------------------------
# sample processing of Azimuth processing for reference labelling
# this part produce the file needed for the actual label transfer

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
# options(Seurat.object.assay.version = "v5")

# read in the dataset -----------------------------------------------------
# this is the query dataset
data.combined <- readRDS(file = "out/object/analysis_R45/04_sobj_filtered_harmony_afterQC.rds")
DimPlot(data.combined,reduction = "umap",group.by = "RNA_snn_res.0.4")

# generate the diet object and save it
# save the coordiantes of the UMAP
data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode") %>%
  write_tsv("out/table/analysis_R45/05_sobj_filtered_harmony_afterQC_coordUMAP.tsv")

DefaultAssay(data.combined) <- "RNA"
object_diet <- DietSeurat(object = data.combined, assays = "RNA")
saveRDS(object_diet,"out/object/analysis_R45/05_sobj_filtered_harmony_afterQC_diet.rds")

# available datasets ------------------------------------------------------
# pick a reference
# load the reference dataset
reference <- LoadReference(path = "data/reference/azimuth/ref_snRNAseq_human_cerebellum")
DimPlot(reference$plot,group.by = "annotation.l1")

# save the reference umap
df_point <- reference$map@reductions$refUMAP@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")
write_tsv(df_point,"out/table/analysis_R45/azimuth_cerebellum_CoordUMAP.tsv")

# save the meta of the ref
df_meta <- reference$map@meta.data %>%
  data.frame() %>%
  rownames_to_column("barcode")
write_tsv(df_meta,"out/table/analysis_R45/azimuth_cerebellum_Metadata.tsv")
