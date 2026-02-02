# AIM ---------------------------------------------------------------------
# Try to repurpose some real data to be used for testing

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(harmony)
library(AnnotationDbi)
library(AnnotationHub)

# read in the object ------------------------------------------------------
test <- readRDS("../../data/c18b60ea-7dbc-4705-a3dc-e29da4e43c68.rds")
test@meta.data
DimPlot(test,split.by = "cell_type",label = T,ncol=4)

# wrangling ---------------------------------------------------------------
# define the assay to pull
table(test$assay)

# define the samples 
as.character(unique(test@meta.data$cell_type))

# subset only the cells of interest
test@meta.data$test_cellid <- test@meta.data %>%
  mutate(test = case_when(cell_type %in% c("erythrocyte") ~ "erythrocyte",
          cell_type %in% c("CD4-positive, alpha-beta T cell",
                           "CD4-positive, alpha-beta memory T cell",
                           "naive thymus-derived CD4-positive, alpha-beta T cell") ~ "CD4 T cell",
          cell_type %in% c("CD8-positive, alpha-beta cytokine secreting effector T cell",
                           "CD8-positive, alpha-beta T cell") ~ "CD8 T cell",
          cell_type %in% c("classical monocyte",
                           "monocyte") ~ "monocyte",
          cell_type %in% c("memory B cell",
                           "naive B cell") ~ "B cells",
          cell_type %in% c("neutrophil") ~ "neutrophil",
          cell_type %in% c("plasma cell") ~ "plasma cell",
          cell_type %in% c("platelet") ~ "platelet",
          cell_type %in% c("mature NK T cell",
                           "type I NK T cell") ~ "NK T cell",
          T~"remove")) %>%
  pull(test)

# labels not included in the final object
# "macrophage"
# "plasmacytoid dendritic cell"
# "plasmablast"                                                
# "non-classical monocyte"
# "hematopoietic stem cell"
# "common myeloid progenitor"
# "granulocyte"
# "basophil"
# "T cell"
# "CD141-positive myeloid dendritic cell"

# perform the filtering over the labels and the assay
test_filter <- subset(test,subset = test_cellid != "remove" & assay == "10x 3' v3")

# add the fake technical metadata
set.seed(21)
donor_pool <- paste0("don_",str_pad(1:6,pad = "0",width = 2))
test_filter$test_donor <- sample(donor_pool,size = dim(test_filter)[2],replace = T)

# Add the corresponding disease condition to be tested
LUT_donor <- data.frame(test_donor = donor_pool) %>%
  mutate(test_disease = case_when(test_donor %in% c("don_01","don_02","don_03")~"ctrl",
                                  T~"dis"))

# add it to the metadata
test_filter$test_disease <- test_filter@meta.data %>%
  left_join(LUT_donor,by="test_donor") %>%
  pull(test_disease)

# check the metadata
table(test_filter$test_cellid,test_filter$test_donor)
table(test_filter$test_cellid,test_filter$test_disease)

# generate the two objects one intact and one with the simulated removal of the neutrophils from disease condition only
# select the cells to be removed
meta_filter <- test_filter@meta.data %>%
  rownames_to_column() %>%
  mutate("barcodes" = rowname)

# randomly select some cells to be remove
set.seed(21)
id_barcodes_remove <- meta_filter %>%
  filter(test_disease == "dis",
         test_cellid == "neutrophil") %>% 
  slice_sample(prop = 0.5) %>%
  pull("barcodes")

meta_filter_full <- meta_filter %>%
  mutate(test_remove = case_when(barcodes %in% id_barcodes_remove ~ "remove",
                                 T ~ "keep"))

# add the metadata for the original object
test_filter$barcodes <- meta_filter_full$barcodes
test_filter$test_remove <- meta_filter_full$test_remove

# save the original object withour remove the cells
test_filter_full <- test_filter

# remove the cells and create a new object
test_filter_remove <- subset(test_filter,subset = test_remove == "keep")

# check the two dataset
table(test_filter_full$test_cellid,test_filter_full$test_donor)
table(test_filter_full$test_cellid,test_filter_full$test_disease)

table(test_filter_remove$test_cellid,test_filter_remove$test_donor)
table(test_filter_remove$test_cellid,test_filter_remove$test_disease)

# extract the expression matrices and rename the genes
LUT_genes <- data.frame(gene_id = rownames(test_filter_full@assays$RNA$counts)) %>%
  left_join(test@assays$RNA@meta.features %>%
              rownames_to_column("gene_id"))

# confirm the source of the genes is the same for both objects, so I can use the same LUT
sum(!(rownames(test_filter_full@assays$RNA$counts) == rownames(test_filter_remove@assays$RNA$counts)))

# change the gene name for both matrices of the version of the data
mat_filter_full <- test_filter_full@assays$RNA$counts
rownames(mat_filter_full) <- LUT_genes$feature_name

mat_filter_remove <- test_filter_remove@assays$RNA$counts
rownames(mat_filter_remove) <- LUT_genes$feature_name

# produce the seurat objects ----------------------------------------------
# notice it is critacal that all matrices have the same dimension
# I need to create a single object to add the cell cycle scoring and other metadata. I decided to keep all the cells and all the genes for now.
sobj_total_full <- CreateSeuratObject(counts = mat_filter_full,
                                      project = "blood_full",
                                      meta.data = test_filter_full@meta.data,
                                      min.cells = 0, min.features = 0) %>%
  # this is needed as the cell cycle scoring is done on the data slot, which would be empty
  Seurat::NormalizeData(verbose = T)

sobj_total_remove <- CreateSeuratObject(counts = mat_filter_remove,
                                        project = "blood_full",
                                        meta.data = test_filter_remove@meta.data,
                                        min.cells = 0, min.features = 0) %>%
  # this is needed as the cell cycle scoring is done on the data slot, which would be empty
  Seurat::NormalizeData(verbose = T)

# add the cell cycle analysis
DefaultAssay(sobj_total_full) <- "RNA"
DefaultAssay(sobj_total_remove) <- "RNA"

# once updated pick the annotation of interest
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sobj_total_full <- CellCycleScoring(sobj_total_full, s.features = s.genes, g2m.features = g2m.genes)
sobj_total_remove <- CellCycleScoring(sobj_total_remove, s.features = s.genes, g2m.features = g2m.genes)

sobj_total_full$percent.mt <- PercentageFeatureSet(sobj_total_full, pattern = "^MT-")
sobj_total_remove$percent.mt <- PercentageFeatureSet(sobj_total_remove, pattern = "^MT-")

sobj_total_full$percent.ribo <- PercentageFeatureSet(sobj_total_full, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
sobj_total_remove$percent.ribo <- PercentageFeatureSet(sobj_total_remove, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")

# add also the percentage of globin. in this dataset it is not meaningful as there is no blood
sobj_total_full$percent.globin <- Seurat::PercentageFeatureSet(sobj_total_full,pattern = "^HB[^(P)]")
sobj_total_remove$percent.globin <- Seurat::PercentageFeatureSet(sobj_total_remove,pattern = "^HB[^(P)]")

# integration processing --------------------------------------------------
# check the scale matrix
sobj_total_full@assays$RNA
sobj_total_remove@assays$RNA
# pull all the genes to scale
# all.genes <- rownames(sobj_total)

# rescale the data for regressing out the sources of variation do not scale all the genes. if needed I can scale them before the heatmap call. for speeding up the computation I will keep 
sobj_total_full <- sobj_total_full %>%
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
  # FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  FindClusters(resolution = 0.4) %>%
  identity()

sobj_total_remove <- sobj_total_remove %>%
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
  # FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  FindClusters(resolution = 0.4) %>%
  identity()

# check the status of dataset preintegration
DimPlot(sobj_total_full,group.by = "test_donor",raster = T)
DimPlot(sobj_total_full,group.by = "donor_id",raster = T)
DimPlot(sobj_total_full,group.by = "test_cellid",raster = T,split.by = "donor_id")

DimPlot(sobj_total_full,group.by = "test_cellid",raster = T,split.by = "test_donor")
DimPlot(sobj_total_full,group.by = "test_cellid",raster = T,label = T)

DimPlot(sobj_total_remove,group.by = "test_donor",raster = T)

# Run Harmony -------------------------------------------------------------
# The simplest way to run Harmony is to pass the Seurat object and specify which variable(s) to integrate out. RunHarmony returns a Seurat object, updated with the corrected Harmony coordinates. Let's set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round.
sobj_total_full_h <- sobj_total_full %>%
  RunHarmony("donor_id", plot_convergence = TRUE)

sobj_total_remove_h <- sobj_total_remove %>%
  RunHarmony("donor_id", plot_convergence = TRUE)

# Downstream analysis -----------------------------------------------------
# Many downstream analyses are performed on low dimensional embeddings, not gene expression. To use the corrected Harmony embeddings rather than PCs, set reduction = 'harmony'. For example, let's perform the UMAP and Nearest Neighbor analyses using the Harmony embeddings.
sobj_total_full_h <- sobj_total_full_h %>%
  RunUMAP(reduction = "harmony", dims = 1:30,return.model = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  # FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  FindClusters(resolution = 0.4) %>%
  identity()

sobj_total_remove_h <- sobj_total_remove_h %>%
  RunUMAP(reduction = "harmony", dims = 1:30,return.model = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  # FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  FindClusters(resolution = 0.4) %>%
  identity()

DimPlot(sobj_total_full_h,group.by = "test_donor",raster = T)
DimPlot(sobj_total_full_h,group.by = "donor_id",raster = T)
DimPlot(sobj_total_full_h,group.by = "test_cellid",raster = T,split.by = "donor_id")

DimPlot(sobj_total_full_h,group.by = "test_cellid",raster = T,split.by = "test_donor")
DimPlot(sobj_total_full_h,group.by = "test_cellid",raster = T,label = T)

DimPlot(sobj_total_remove_h,group.by = "test_donor",raster = T)
DimPlot(sobj_total_remove_h,group.by = "test_cellid",raster = T,label = T)

# save objects ------------------------------------------------------------
saveRDS(sobj_total_full_h,"../../out/object/sobj_total_full_h.rds")
saveRDS(sobj_total_remove_h,"../../out/object/sobj_total_remove_h.rds")
