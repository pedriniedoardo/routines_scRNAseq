# AIM ---------------------------------------------------------------------
# read in the reference dataset from cellxgene
# produce an object usable as reference for the annotaion of the query dataset

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(schard)
library(Azimuth)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
# options(Seurat.object.assay.version = "v3")

# in cause of big references increase the amount of RAM per worker
options(future.globals.maxSize = 15 * 1000 * 1024^2) # 15 GB
# print the default number of dimensions for Azimuth
options()$Azimuth.map.ndims

# read in the data --------------------------------------------------------
scobj <- schard::h5ad2seurat('data/reference/1201b573-417d-4a9d-86f7-c43cc2bc36e3.h5ad')

# explore the object
scobj@meta.data

# filter only the barcodes of interest
# only from control samples
# only from fully developed samples
# only from known cell types
scobj@meta.data$disease %>% table()
scobj@meta.data$development_stage %>% table()
scobj@meta.data$cell_type %>% table()

scobj@meta.data %>%
  filter(cell_type %in% c("Purkinje cell")) %>%
  group_by(development_stage) %>%
  summarise(n = n())

# apply the filter
scobj_filter <- scobj %>%
  subset(subset = disease == "normal") %>%
  subset(subset = development_stage %in% c("42-year-old stage","44-year-old stage","46-year-old stage","52-year-old stage")) %>%
  subset(subset = cell_type != "unknown")

# regenerate the object
DimPlot(scobj_filter,group.by = "cell_type")

# generate the azimuth reference object -----------------------------------

# processing --------------------------------------------------------------
# define the ident for the object
Idents(object = scobj_filter) <- "cell_type"

# the object needs to have a model saved, also it needs to have the sct transformed data
ref2 <- SCTransform(scobj_filter, method = "glmGamPoi", verbose = T) %>%
  # it is recommended to have 50 dimensions for a reliable map transfering: options()$Azimuth.map.ndims
  RunPCA(verbose = T,npcs = 50) %>%
  RunUMAP(dims = 1:30, return.model = TRUE)

# after SCT the dimensionality reduction is recomputed
DimPlot(ref2,group.by = "cell_type",label = T)
DimPlot(ref2,group.by = "development_stage",label = T)
DimPlot(ref2,group.by = "cell_type",split.by = "development_stage",label = T)

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
SaveAnnoyIndex(object = ref_final[["refdr.annoy.neighbors"]], file = file.path("data/reference/azimuth/ref_snRNAseq_human_cerebellum/","idx.annoy"))
saveRDS(object = ref_final, file = file.path("data/reference/azimuth/ref_snRNAseq_human_cerebellum/", "ref.Rds"))
saveRDS(object = full.ref, file = file.path("data/reference/azimuth/ref_snRNAseq_human_cerebellum/", "fullref.Rds"))

# test reading azimuth reference ------------------------------------------
# try to load the new reference object created

# if needed reduce the dimensions of the object before loading it
# options(Azimuth.map.ndims = 30)
reference <- LoadReference(path = "data/reference/azimuth/ref_snRNAseq_human_cerebellum/")
DimPlot(reference$plot,group.by = "annotation.l1",label = T)

