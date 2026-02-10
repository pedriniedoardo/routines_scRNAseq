# AIM ---------------------------------------------------------------------
# try to run harmony by merging the matrices from the individula objects. this will allow the skipping of the regular integration. the regular integration is needed as to run harmony the matrices should have the same number/order of the genes.
# for this step I have used the 500 5000 10 version of the filtering

# renv integration --------------------------------------------------------

# to load the packages
source(".Rprofile")

# in the config specify the following
# renv_library_path: "renv/library/linux-rocky-9.5/R-4.5/x86_64-conda-linux-gnu"

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
library(homologene)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
# load the LUT
# LUT <- read_csv("../../data/LUT_samples.csv")

# read in the list of objects. use the filtered dataset for the singlets only
data.list <- readRDS("out/object/analysis_R45/03_list_sobj_filtered_afterQC.rds")
label <- "afterQC"

# read in the mouse - human annotation
# homologeneData2_250121 <- read_tsv("../../data/Homologene_250121") %>%
#   as.data.frame()

# merge the reads in a single table to allow harmony processing -----------
# merge the individual objetcs (already filtered as individula objects) to create a single count matrix
data.list.id <- str_remove_all(names(data.list),pattern = ".rds")
data.combined.all <- merge(data.list[[1]], y = data.list[-1], add.cell.ids = data.list.id, project = "test_V5")

# check the size of the dataset
data.combined.all

# confirm the total number of cells
lapply(data.list,function(x){
  dim(x@assays$RNA$counts)[2]
}) %>%
  unlist() %>%
  sum()

# join the layers to be able to create a single layer
# this step is needed for the new assay 5 structure
# see the number of layers before joining
names(data.combined.all@assays$RNA@layers)
data.combined.all[["RNA"]] <- JoinLayers(data.combined.all[["RNA"]])

# see the number of layers after joining
names(data.combined.all@assays$RNA@layers)

# generate a minimal metadata
# remove the provious clustering resutls from the individual objects
meta_minimal <- data.combined.all@meta.data %>%
  dplyr::select(-c(seurat_clusters,contains("RNA_snn_res")))

# create the seurat object ------------------------------------------------
# notice it is critical that all matrices have the same dimension
# I need to create a single object to add the cell cycle scoring and other metadata. I decided to trimm further the dataset for genes content
# in case you want to keep all the genes from the merged object do not apply any filter for min.cells and min.features
# NOTICE!!! using data.combined.all@assays$RNA$counts retains the gene names in the matrix
sobj_total <- CreateSeuratObject(counts = data.combined.all@assays$RNA$counts,
                                 project = "Maltecca_2552_V5",
                                 meta.data = meta_minimal,
                                 # remove the low expressing genes
                                 min.cells = 20,
                                 # keep all the cells in this case, the filtering has already been performed
                                 min.features = 0) %>%
  # this is needed as the cell cycle scoring is done on the data slot, which would be empty
  Seurat::NormalizeData(verbose = T)

# after creating the object I do not need the list and the merged matrix, free up some space
remove(data.list,data.combined.all)
gc()

# add the cell cycle analysis
DefaultAssay(sobj_total) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# s.genes.mouse <- homologene(s.genes, inTax = 9606, outTax = 10090,db = homologeneData2_250121) %>%
#   pull(`10090`) %>%
#   unique()
# g2m.genes.mouse <- homologene(g2m.genes, inTax = 9606, outTax = 10090,db = homologeneData2_250121) %>%
#   pull(`10090`) %>%
#   unique()

# sobj_total <- CellCycleScoring(sobj_total, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse)
sobj_total <- CellCycleScoring(sobj_total, s.features = s.genes, g2m.features = g2m.genes)

# in this case the covariates are alredy present in the metadata so I can skip the routine below for calcularing the QC covariates

# # add also the percentage of globin. in this dataset it is not meaningful as there is no blood
# # datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
# sobj_total$percent.mt.integration <- PercentageFeatureSet(sobj_total, pattern = "^mt-")
# sobj_total$percent.ribo.integration <- PercentageFeatureSet(sobj_total, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")
# # add also the percentage of globin. in this dataset it is not meaningful as there is no blood
# sobj_total$percent.globin.integration <- Seurat::PercentageFeatureSet(sobj_total,pattern = "^Hb[^(p)]")
# 
# # confirm the proportions are the same values
# all.equal(sobj_total$percent.mt.integration,sobj_total$percent.mt)
# all.equal(sobj_total$percent.ribo.integration,sobj_total$percent.ribo)
# all.equal(sobj_total$percent.globin.integration,sobj_total$percent.globin)

# standard processing -----------------------------------------------------
# check the scale matrix
sobj_total@assays$RNA$scale.data

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
# sobj_total_h@assays$RNA@counts[1:20,1:10]
# sobj_total_h@assays$RNA@data[1:20,1:10]
# sobj_total_h@assays$RNA@scale.data[1:20,1:10]
# 
# dim(sobj_total_h@assays$RNA@counts)
# dim(sobj_total_h@assays$RNA@data)
# dim(sobj_total_h@assays$RNA@scale.data)
sobj_total_h@assays$RNA$counts[1:20,1:10]
sobj_total_h@assays$RNA$data[1:20,1:10]
sobj_total_h@assays$RNA$scale.data[1:20,1:10]

dim(sobj_total_h@assays$RNA$counts)
dim(sobj_total_h@assays$RNA$data)
dim(sobj_total_h@assays$RNA$scale.data)

DimPlot(sobj_total_h,group.by = "orig.ident",raster = F)

# save the object
saveRDS(sobj_total_h,paste0("out/object/analysis_R45/04_sobj_filtered_harmony_",label,".rds"))
