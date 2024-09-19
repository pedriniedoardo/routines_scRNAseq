# AIM ---------------------------------------------------------------------
# try to run harmony by merging the matrices from the individula objects. this will allow the skipping of the regular integration. the regular integration is needed as to run harmony the matrices should have the same number/order of the genes.
# for this step I decided to use the 01000_6000_15 version of the filtering. the cells of intererest are well mainteined

# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)
library(harmony)
library(ggrepel)
library(ComplexHeatmap)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v3")

# read in the data --------------------------------------------------------
# load the LUT
# LUT <- read_csv("../../data/LUT_samples.csv")

# read in the list of objects. use the filtered dataset for the singlets only
data.list <- readRDS("../out/test_introns/object/list_datasc_fix_filter_norm_doubletSinglet_SoupX_01000_06000_15_V3.rds")

# merge the reads in a single table to allow harmony processing -----------
# merge the individual objetcs (already filtered as individula objects) to create a single count matrix
data.list.id <- str_remove_all(names(data.list),pattern = ".rds")
data.combined.all <- merge(data.list[[1]], y = data.list[-1], add.cell.ids = data.list.id, project = "test_V5")

# check the size of the dataset
data.combined.all

# old AssayV3
# confirm the total number of cells
lapply(data.list,function(x){
  dim(x@assays$RNA@counts)[2]
}) %>%
  unlist() %>%
  sum()

# create the seurat object ------------------------------------------------
# notice it is critical that all matrices have the same dimension
# I need to create a single object to add the cell cycle scoring and other metadata. I decided to trimm further the dataset for genes content
# in case you want to keep all the genes from the merged object do not apply any filter for min.cells and min.features

# old AssayV3
# notice that for the old assay version there is non need to normalize the data before the cell cycle scoring. In this case I am doing it now to keep the most similar workflow.
sobj_total <- CreateSeuratObject(counts = data.combined.all@assays$RNA@counts,
                                 project = "test_V5",
                                 meta.data = data.combined.all@meta.data,
                                 min.cells = 20, min.features = 200) %>%
  # NOTICE: If I normalize before the cell cycle assignament I obtain a different call for the cell cycle identity.
  Seurat::NormalizeData(verbose = T)

# after creating the object I do not need the list and the merged matrix, free up some space
remove(data.list,data.combined.all)
gc()

# add the cell cycle analysis
DefaultAssay(sobj_total) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sobj_total <- CellCycleScoring(sobj_total, s.features = s.genes, g2m.features = g2m.genes)
sobj_total$percent.mt.integration <- PercentageFeatureSet(sobj_total, pattern = "^MT-")
sobj_total$percent.ribo.integration <- PercentageFeatureSet(sobj_total, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
# add also the percentage of globin. in this dataset it is not meaningful as there is no blood
sobj_total$percent.globin.integration <- Seurat::PercentageFeatureSet(sobj_total,pattern = "^Hb[^(p)]")

# if needed integrate some new metadata at this point ---------------------
# # add all the original metadata 
# meta_new <- sobj_total@meta.data %>% 
#   rownames_to_column("barcode")
# 
# meta_new_total <- left_join(meta_new,df_meta %>% dplyr::select(barcode,seurat_clusters_ref),"barcode",suffix=c(".harmony",".cca")) %>% 
#   left_join(LUT,by = c("orig.ident"="official_id")) %>% 
#   column_to_rownames("barcode")
# 
# # update the meta
# sobj_total@meta.data <- meta_new_total

# standard processing -----------------------------------------------------
# check the scale matrix
sobj_total@assays$RNA$scale.data

# old AssayV3
# sobj_total@assays$RNA@scale.data

# I generally do not scale all the features in the object to reduce computation time and size of the object.
# pull all the genes to scale
# all.genes <- rownames(sobj_total)

# rescale the data for regressing out the sources of variation do not scale all the genes.
# if needed for some plots, I can scale them before the heatmap call. for speeding up the computation I will scale only the HVF
sobj_total <- sobj_total %>%
  # skip the normalizatio that has been already performed at the beginning
  # Seurat::NormalizeData(verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>%
  # I can scale the missing features afterwards now focus on the highly variable one for speed purposes
  ScaleData(vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"), verbose = T) %>% 
  # run this if you want to scale all the variables
  # ScaleData(vars.to.regress = c("percent.mt.harmony","nCount_RNA.harmony","S.Score.harmony","G2M.Score.harmony"), verbose = T,features = all.genes) %>% 
  RunPCA(npcs = 30, verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:30,return.model = TRUE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  identity()

# check the status of dataset preintegration
DimPlot(sobj_total,group.by = "orig.ident",raster = T)

# Run Harmony -------------------------------------------------------------
# The simplest way to run Harmony is to pass the Seurat object and specify which variable(s) to integrate out. RunHarmony returns a Seurat object, updated with the corrected Harmony coordinates. Let's set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round.
sobj_total_h <- sobj_total %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

# Harmony with two or more covariates
# Do the same with your Seurat object:
# seuratObject <- RunHarmony(seuratObject, c("dataset", "donor", "batch_id"))
# To directly access the new Harmony embeddings, use the Embeddings command.
harmony_embeddings <- Embeddings(sobj_total_h, 'harmony')
harmony_embeddings[1:5, 1:5]
# Let's make sure that the datasets are well integrated in the first 2 dimensions after Harmony.
# DimPlot(object = sobj_total_h, reduction = "harmony", pt.size = .1, group.by = "sample_id")

# Downstream analysis -----------------------------------------------------
# Many downstream analyses are performed on low dimensional embeddings, not gene expression. To use the corrected Harmony embeddings rather than PCs, set reduction = 'harmony'. For example, let's perform the UMAP and Nearest Neighbor analyses using the Harmony embeddings.
sobj_total_h <- sobj_total_h %>%
  RunUMAP(reduction = "harmony", dims = 1:30,return.model = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  # FindClusters(resolution = 0.5) %>%
  identity()

# verify that all the relevant slots are filled
sobj_total_h@assays$RNA@counts[1:20,1:10]
sobj_total_h@assays$RNA@data[1:20,1:10]
sobj_total_h@assays$RNA@scale.data[1:20,1:10]

dim(sobj_total_h@assays$RNA@counts)
dim(sobj_total_h@assays$RNA@data)
dim(sobj_total_h@assays$RNA@scale.data)

DimPlot(sobj_total_h,group.by = "ID",raster = F)

# save the object
saveRDS(sobj_total_h,"../out/test_introns/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15_V3.rds")
