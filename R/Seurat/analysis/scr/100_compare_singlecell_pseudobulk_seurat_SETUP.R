# AIM ---------------------------------------------------------------------
# the aim of this test is to compare seurat's implementaiton of pseudobulk analysis, to the default process of DGE
# the reference of the test is presented here: https://satijalab.org/seurat/articles/de_vignette.html
# this first script just process the dataset to add important metadata to the object

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(tidyverse)
library(harmony)

# read in the data --------------------------------------------------------
# check available dataset
# SeuratData::AvailableData()

# install the dataset
# SeuratData::InstallData("ifnb")

# load the dataset
ifnb <- SeuratData::LoadData("ifnb")
table(ifnb@meta.data$orig.ident)

# wrangling ---------------------------------------------------------------
# add the donor ID to the metadata
# download.file("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye1.ctrl.8.10.sm.best",destfile = "../data/misc/ye1.ctrl.8.10.sm.best.tsv")
# download.file("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye2.stim.8.10.sm.best",destfile = "../data/misc/ye2.stim.8.10.sm.best.tsv")
# ctrl <- read.table("../data/misc/ye1.ctrl.8.10.sm.best.tsv",head = T, stringsAsFactors = F)
# stim <- read.table("../data/misc/ye2.stim.8.10.sm.best.tsv",head = T, stringsAsFactors = F)
ctrl <- read_tsv("../data/misc/ye1.ctrl.8.10.sm.best.tsv")
stim <- read_tsv("../data/misc/ye2.stim.8.10.sm.best.tsv")

# merge the two tables
info <- bind_rows(ctrl, stim)
df_info <- info %>%
  # rename the cell IDs by substituting the '-' into '.'
  mutate(BARCODE = str_replace_all(BARCODE,pattern = "\\-", replacement = "\\.")) %>%
  # only keep the cells with high-confidence sample ID
  filter(str_detect(BEST,"SNG"))

# identify the duplicated barcodes
barcodes_remove <- df_info %>%
  mutate(dup = duplicated(BARCODE)) %>%
  filter(dup == T) %>%
  pull(BARCODE)

# remove all the barcodes that are duplicated
df_info_final <- df_info %>%
  filter(!(BARCODE %in% barcodes_remove)) %>%
  column_to_rownames("BARCODE") %>%
  select(BEST) %>%
  dplyr::rename(donor_id = BEST)

# add the meta to the object
ifnb2 <- AddMetaData(ifnb, metadata = df_info_final)

# confirm the barcodes id are matching
# ifnb2@meta.data %>%
#   rownames_to_column("barcode") %>%
#   left_join(df_info_final %>% rownames_to_column("barcode"),by = "barcode") %>%
#   mutate(test = donor_id.x == donor_id.y) %>%
#   filter(!is.na(test)) %>%
#   filter(test == F)

# give the donor unknown label to the cell without a donor id
vec_donor_id <- ifnb2@meta.data %>%
  mutate(donor_id_fix = case_when(is.na(donor_id) ~ "unknown",
                                  T ~ donor_id)) %>%
  pull(donor_id_fix)
# add the new metadata to the object
ifnb2$donor_id_fix <- vec_donor_id

# add the donor+stimulation variable
ifnb2$stim_donor <- paste0(ifnb2$stim,"_",ifnb2$donor_id_fix)

# filter the dataset
ifnb_final <- subset(ifnb2,subset = donor_id_fix != "unknown")

# notice that in this dataset the sample is parid with both stimulated and unstimulated donors
ifnb_final@meta.data %>%
  group_by(orig.ident,donor_id_fix) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = orig.ident,values_from = n,values_fill = 0)

# how many cells are in the dataset per annotation
ifnb_final@meta.data %>%
  group_by(seurat_annotations) %>%
  summarise(n = n())

ifnb_final@meta.data %>%
  group_by(seurat_annotations, stim) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = stim,values_from = n,values_fill = 0)

ifnb_final@meta.data %>%
  group_by(seurat_annotations, stim,donor_id) %>%
  summarise(n = n()) %>%
  filter(seurat_annotations %in% c("CD14 Mono")) %>%
  pivot_wider(names_from = stim,values_from = n,values_fill = 0)

# standard pre-processing of the sample -----------------------------------
sobj_total <- ifnb_final %>%
  # this is needed as the cell cycle scoring is done on the data slot, which would be empty
  Seurat::NormalizeData(verbose = T)

# add the cell cycle analysis
DefaultAssay(sobj_total) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# add cell cycle scoring
sobj_total <- CellCycleScoring(sobj_total, s.features = s.genes, g2m.features = g2m.genes)
# add prop of mito reads, but it seems the mito genes have been hard filtered out of the matrix
sobj_total$percent.mt <- PercentageFeatureSet(sobj_total, pattern = "^MT-")
sobj_total$percent.ribo <- PercentageFeatureSet(sobj_total, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
# add also the percentage of globin. in this dataset it is not meaningful as there is no blood
sobj_total$percent.globin <- Seurat::PercentageFeatureSet(sobj_total,pattern = "^HB[^(P)]")

# check the scale matrix
sobj_total@assays$RNA$scale.data

# 
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
DimPlot(sobj_total,group.by = "seurat_annotations",raster = T)
DimPlot(sobj_total,group.by = "donor_id_fix",raster = T)

# Run Harmony -------------------------------------------------------------
# considering the stricture of teh dataset try to regress out the stimulation + donor
sobj_total_h <- sobj_total %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

sobj_total_h2 <- sobj_total %>%
  RunHarmony("stim_donor", plot_convergence = TRUE)

# Downstream analysis -----------------------------------------------------
# Many downstream analyses are performed on low dimensional embeddings, not gene expression. To use the corrected Harmony embeddings rather than PCs, set reduction = 'harmony'. For example, let's perform the UMAP and Nearest Neighbor analyses using the Harmony embeddings.
sobj_total_h <- sobj_total_h %>%
  RunUMAP(reduction = "harmony", dims = 1:30,return.model = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  # FindClusters(resolution = 0.5) %>%
  identity()

sobj_total_h2 <- sobj_total_h2 %>%
  RunUMAP(reduction = "harmony", dims = 1:30,return.model = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  # FindClusters(resolution = 0.5) %>%
  identity()

(DimPlot(sobj_total_h,group.by = "orig.ident",raster = T)+ggtitle("int sample")) + (DimPlot(sobj_total_h2,group.by = "orig.ident",raster = T)+ggtitle("int sample_stim"))
DimPlot(sobj_total_h,group.by = "seurat_annotations",raster = T) + DimPlot(sobj_total_h2,group.by = "seurat_annotations",raster = T)
DimPlot(sobj_total_h,group.by = "donor_id_fix",raster = T) + DimPlot(sobj_total_h2,group.by = "donor_id_fix",raster = T)

DimPlot(sobj_total_h,split.by = "donor_id_fix",raster = T)
DimPlot(sobj_total_h,group.by = "seurat_annotations",raster = T,label = T)

# save the object
saveRDS(sobj_total_h,"../out/object/100_ifnb_stim.rds")
saveRDS(sobj_total_h2,"../out/object/100_ifnb_DonorStim.rds")

