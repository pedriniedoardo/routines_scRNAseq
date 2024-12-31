# AIM ---------------------------------------------------------------------
# this vignette will show how to use the SCType package to annotate the cell types in the dataset

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ggraph)
library(igraph)
library(data.tree)
library(HGNChelper)
library(openxlsx)

# functions ---------------------------------------------------------------
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R",destfile = "scr/SCType/gene_sets_prepare.R")
source("scr/SCType/gene_sets_prepare.R")
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R",destfile = "scr/SCType/sctype_score_.R")
source("scr/SCType/sctype_score_.R")
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R",destfile = "scr/SCType/auto_detect_tissue_type.R")
source("scr/SCType/auto_detect_tissue_type.R")

# download the database file
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",destfile = "../data/ScTypeDB_full.xlsx")
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",destfile = "../data/ScTypeDB_short.xlsx")

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v3")


# read in the dataset -----------------------------------------------------
scobj <- readRDS("../out/test_introns/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15_V3.rds")
# define which is the resolution to work with
# plot the UMAP with all the resolutions runs
id_resolution <- str_subset(colnames(scobj@meta.data),pattern = "RNA_snn_res") %>%
  str_subset(negate = T,pattern = ".paper") %>%
  sort()

list_plot <- lapply(id_resolution,function(x){
  plot <- DimPlot(scobj,
                  reduction = "umap",
                  group.by = x,
                  label = T,
                  raster = T)
  return(plot)
})

wrap_plots(list_plot)
ggsave("../out/test_introns/plot/03_UMAPCluster_resolutions_V3.pdf",width = 20,height = 15)

# pick the resolution
ident_resolution <- "RNA_snn_res.0.2"

# rename the seurat cluster varaible with the picked resolution
scobj$seurat_clusters <- scobj@meta.data[[ident_resolution]]
Idents(scobj) <- "seurat_clusters"

# confirm the identity of the object
DimPlot(scobj,label = T,raster = T)

# Confirm the object has a scaled matrix in it. this is the one that will bu used for the scoring of the cell types
# this means that if a gene is not in the scaled matrix it will not contribute to the label scoring
# data.combined@assays$RNA@scale.data[1:10,1:10]
# dim(data.combined@assays$RNA@scale.data)
scobj@assays$RNA@scale.data[1:10,1:10]
dim(scobj@assays$RNA@scale.data)

# read in the database file
db_ <- "../data/ScTypeDB_full.xlsx"
sample_db <- readxl::read_excel(db_)

# pick the tissue of reference for the anntoaiton
# e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
tissue <- "Brain"

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)
str(gs_list)

# get cell-type by cell matrix score
es.max <- sctype_score(scRNAseqData = scobj[["RNA"]]@scale.data,
                       scaled = TRUE, 
                       gs = gs_list$gs_positive,
                       gs2 = gs_list$gs_negative) 

es.max[1:10,1:10]

# wrangling ---------------------------------------------------------------
# after running the annotation in theory all the cells can have a specific annotation
df_score <- es.max %>% 
  data.frame() %>% 
  rownames_to_column("annotation") %>% 
  pivot_longer(names_to = "barcode",values_to = "score",-annotation) %>% 
  mutate(barcode = str_replace(barcode,pattern = "\\.",replacement = "-"))

df_score

# add the new annotation to the original object and show the umap
df_meta <- scobj@meta.data %>% 
  rownames_to_column("barcode")

head(df_meta)

# for each cluster in the dataset sum the scores per cell for each cluster
# cl <- 0
df_score_cluster <- lapply(unique(Idents(scobj)),function(cl){
  # pull the cells belonging to the cluster
  id <- df_meta %>%
    
    filter(seurat_clusters %in% cl) %>% 
    pull(barcode)
  
  # sum the score from eah cell for each category
  df_tot <- df_score %>% 
    filter(barcode %in% id) %>% 
    # add the cluster information
    mutate(cluster = cl) %>% 
    # add the total number of cells
    # mutate(n_cell = length(id)) %>% 
    group_by(annotation,cluster) %>% 
    summarise(sum_score = sum(score),
              n_cell = n())
  
  return(df_tot)
  
}) %>% 
  bind_rows()

df_score_cluster

# pull the top score per cluster
df_score_cluster_top <- df_score_cluster %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = sum_score) %>% 
  # set low-confident (low ScType score) clusters to "unknown"
  mutate(annotation_confident = case_when(sum_score < n_cell/4 ~"Unknown",
                                          T~annotation)) %>% 
  ungroup()

df_score_cluster_top

# add the annotation to the original dataset
df_meta_full <- left_join(df_meta,df_score_cluster_top,c("seurat_clusters"="cluster"))
head(df_meta_full)

# add it back to the object
scobj@meta.data <- df_meta_full %>% 
  column_to_rownames("barcode")

head(scobj@meta.data)

# plotting ----------------------------------------------------------------
DimPlot(scobj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'annotation',raster = T)
ggsave("../out/test_introns/plot/03_UMAP_SCType_V3.pdf",width = 10,height = 7)

# save the annotated object -----------------------------------------------
saveRDS(scobj,"../out/test_introns/object/03_SCtypeAnnotation_V3.rds")
