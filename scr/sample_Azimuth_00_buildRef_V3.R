# AIM ---------------------------------------------------------------------
# sample routine for the generation of the reference object for Azimuth

# libraries ---------------------------------------------------------------
# library(Matrix)
library(Seurat)
library(Azimuth)
library(presto)
library(tidyverse)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v3")
# in cause of big references increase the amount of RAM per worker
# options(future.globals.maxSize = 15 * 1000 * 1024^2) # 15 GB

# processing --------------------------------------------------------------
# args <- commandArgs(trailingOnly = TRUE)
# read in the object for which I want to build a reference
ref <- readRDS("../data/data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")
DimPlot(ref,group.by = "seurat_clusters",label = T)

# define the ident for the object
Idents(object = ref) <- "seurat_clusters"

# the object needs to have a model saved, also it needs to have the sct transformed data
ref_SCT <- SCTransform(ref, method = "glmGamPoi", verbose = T)
ref2 <- RunUMAP(object = ref_SCT,dims = 1:10, return.model = TRUE)

# after SCT the dimensionality reduction is recomputed
DimPlot(ref2,group.by = "seurat_clusters",label = T)

# routine to save the reference for azimuth -------------------------------
full.ref <- ref2

# define the colors for the annotaion
full.ref$annotation.l1 <- Idents(object = full.ref)
colormap <- list(annotation.l1 = CreateColorMap(object = ref2, seed = 2))
colormap[["annotation.l1"]] <- colormap[["annotation.l1"]][sort(x = names(x = colormap[["annotation.l1"]]))]

# show the colors
scales::show_col(colormap$annotation.l1)

# Check the PCA reduction information
print(full.ref[['pca']])
num_computed_dims <- ncol(full.ref[['pca']])
print(paste("Number of computed PC dimensions to set is 1:", num_computed_dims))

# build the object
ref_final <- AzimuthReference(
  object = full.ref,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "SCT",
  metadata = c("annotation.l1"),
  dims = 1:num_computed_dims,
  # dims = 1:50,
  # k.param = 31,
  colormap = colormap,
  reference.version = "1.0.0"
)

DimPlot(ref_final)

# save the reference object ----------------------------------------------
SaveAnnoyIndex(object = ref_final[["refdr.annoy.neighbors"]], file = file.path("../data/azimuth/ref_BSrun01run02/","idx.annoy"))
saveRDS(object = ref_final, file = file.path("../data/azimuth/ref_BSrun01run02/", "ref.Rds"))
saveRDS(object = full.ref, file = file.path("../data/azimuth/ref_BSrun01run02/", "fullref.Rds"))

# test reading azimuth reference ------------------------------------------
# try to load the new reference object created
reference <- LoadReference(path = "../data/azimuth/ref_BSrun01run02/")
DimPlot(reference$plot,group.by = "annotation.l1",label = T)
