# AIM ---------------------------------------------------------------------
# attempt general annotation of the dataset using SCType

# LIBRARIES and FUNCTIONS -------------------------------------------------
library(Seurat)
library(tidyverse)
library(ggraph)
library(igraph)
library(data.tree)
library(HGNChelper)
library(cowplot)
library(ComplexHeatmap)

# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R",destfile = "scr/00_gene_sets_prepare.R")
source("scr/SCType/00_gene_sets_prepare.R")
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R",destfile = "scr/00_sctype_score_.R")
source("scr/SCType/00_sctype_score_.R")
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R",destfile = "scr/00_auto_detect_tissue_type.R")
source("scr/SCType/00_auto_detect_tissue_type.R")

# download the database file
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",destfile = "../../data/ScTypeDB_full.xlsx")

# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",destfile = "../../data/ScTypeDB_short.xlsx")

# read in the data --------------------------------------------------------
# read in the list of single objects
list_datasc_pro <- readRDS("../../../out/R_analysis/R_analysis/object/03_list_datasc_fixed_pro.rds")
DimPlot(list_datasc_pro$test_neuron_token,label = T,group.by = "coarse_cell_type") + 
  DimPlot(list_datasc_pro$test_neuron_auto,label = T,group.by = "coarse_cell_type")

# read in the annotation information
db_ <- "../../data/ScTypeDB_full.xlsx"
tissue <- "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)
str(gs_list)

# annotation --------------------------------------------------------------
# scobj <- list_datasc_pro$test_neuron_auto
# loop the annotation on each individual file
list_datasc_pro_SCType <- lapply(list_datasc_pro,function(scobj){
  id_sample <- unique(scobj@meta.data$orig.ident)
  print(id_sample)
  
  # DimPlot(scobj,label = T)
  
  # scobj@assays$RNA@scale.data[1:10,1:10]
  # dim(scobj@assays$RNA@scale.data)
  
  # assign cell types to each cluster:
  # get cell-type by cell matrix
  es.max <- sctype_score(scRNAseqData = scobj[["RNA"]]$scale.data, scaled = TRUE, 
                         gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  
  es.max[1:10,1:10]
  
  # wrangling ---------------------------------------------------------------
  # after running the annotation all the cells can have a specific annotation
  df_score <- es.max %>% 
    data.frame() %>% 
    rownames_to_column("annotation") %>% 
    pivot_longer(names_to = "barcode",values_to = "score",-annotation) %>% 
    mutate(barcode = str_replace(barcode,pattern = "\\.",replacement = "-"))
  
  # df_score
  # add the new annotation to the original object and show the umap
  df_meta <- scobj@meta.data %>% 
    rownames_to_column("barcode")
  
  # head(df_meta)
  
  # for each cluster in the dataset sum the scores per cell for each cluster
  df_score_cluster <- lapply(unique(df_meta$seurat_clusters),function(cl){
    # pull the cells belonging to the clster
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
  
  # df_score_cluster
  # pull the top score per cluster
  df_score_cluster_top <- df_score_cluster %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = sum_score) %>% 
    # set low-confident (low ScType score) clusters to "unknown"
    mutate(annotation_confident = case_when(sum_score < n_cell/4 ~"Unknown",
                                            T~annotation)) %>% 
    ungroup()
  
  # df_score_cluster_top
  
  # add the annotation to the original dataset
  df_meta_full <- left_join(df_meta,df_score_cluster_top,c("seurat_clusters"="cluster"))
  # head(df_meta_full)
  
  # table summaries ---------------------------------------------------------
  # summarise the number of cells per sample
  a <- df_meta_full %>%
    group_by(orig.ident, annotation) %>%
    tally(name = "Number_of_cells") %>%
    mutate(Percentage = (Number_of_cells/sum(Number_of_cells))*100) 
  
  # save the meta of the annotation and the summaries
  write_tsv(a,paste0("../../out/R_analysis/table/04_Summary_SCType_",id_sample,".tsv"))
  write_tsv(df_meta_full,paste0("../../out/R_analysis/table/04_Meta_SCType_",id_sample,".tsv"))
  
  # add it back to the object
  scobj@meta.data <- df_meta_full %>% 
    column_to_rownames("barcode")
  
  # head(scobj@meta.data)
  
  # plotting ----------------------------------------------------------------
  # DimPlot(scobj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'annotation_confident')
  um <- DimPlot(scobj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'annotation')
  
  # save only the individual plot
  ggsave(plot = um,paste0("../../out/R_analysis/plot/04_SCType_UMAP_",id_sample,".pdf"),width = 12,height = 9)
  ggsave(plot = um,paste0("../../out/R_analysis/plot/04_SCType_UMAP_",id_sample,".png"),width = 12,height = 9)
  
  # save the annotation and the cluster id
  um2 <- DimPlot(scobj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'seurat_clusters')
  plot01 <- um2 | um
  ggsave(plot = plot01,paste0("../../out/R_analysis/plot/04_SCType_cluster_UMAP_",id_sample,".pdf"),width = 20,height = 9)
  ggsave(plot = plot01,paste0("../../out/R_analysis/plot/04_SCType_cluster_UMAP_",id_sample,".png"),width = 20,height = 9)
  
  # plot the proportions ----------------------------------------------------
  # extract the meta after annotation
  meta <- scobj@meta.data %>% 
    rownames_to_column("barcodes")
  
  # save the summary of cluster annotation
  pC <- meta %>%
    group_by(seurat_clusters, orig.ident) %>%
    tally() %>%
    ggplot(., aes(x = orig.ident, y = n, fill=seurat_clusters)) +
    geom_bar(stat = "identity", position="fill") +
    ggtitle("Seurat clusters") +
    theme_cowplot() +
    labs(x = "SampleID", y = "fraction of cells") 
  
  # save only the individual plot
  ggsave(plot = pC,paste0("../../out/R_analysis/plot/04_SCTypecluster_prop_",id_sample,".pdf"),width = 5,height = 9)
  ggsave(plot = pC,paste0("../../out/R_analysis/plot/04_SCTypecluster_prop_",id_sample,".png"),width = 5,height = 9)
  
  # save the summary of SingleR annotation
  pN <- meta %>%
    group_by(annotation, orig.ident) %>%
    tally() %>%
    ggplot(., aes(x = orig.ident, y = n, fill=annotation)) +
    geom_bar(stat = "identity", position="fill") +
    ggtitle("SCType predictions") +
    theme_cowplot() +
    labs(x = "SampleID", y = "fraction of cells") 
  
  # save only the individual plot
  ggsave(plot = pN,paste0("../../out/R_analysis/plot/04_SCType_prop_",id_sample,".pdf"),width = 5,height = 9)
  ggsave(plot = pN,paste0("../../out/R_analysis/plot/04_SCType_prop_",id_sample,".png"),width = 5,height = 9)
  
  # save the summary annotation for cluster and SingleR
  plot02 <- pC | pN
  ggsave(plot = plot02,paste0("../../out/R_analysis/plot/04_SCType_cluster_prop_",id_sample,".pdf"),width = 10,height = 9)
  ggsave(plot = plot02,paste0("../../out/R_analysis/plot/04_SCType_cluster_prop_",id_sample,".png"),width = 10,height = 9)
  
  # save the annotated object -----------------------------------------------
  return(scobj)
})

# save the list of the object after SingleR annotation
saveRDS(list_datasc_pro_SCType,"../../out/R_analysis/object/04_list_datasc_pro_SCType.rds")

# check concordance for the annotation ------------------------------------
# use the jaccard score to check the cross simularity of the annotations
# define the jaccard score function
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

# test
a <- c('potato', 'tomotto', 'chips', 'baloon')
b <- c('car', 'chips', 'bird', 'salt')

jaccard(a, b)

df_meta_full2 <- list_datasc_pro_SCType$test_neuron_auto@meta.data %>%
  rownames_to_column("barcode")

#
# build the dataset for the correlatino plot
df_crossing <- tidyr::crossing(SCType_anno = unique(list_datasc_pro_SCType$test_neuron_auto@meta.data$annotation_confident),
                               cellranger_anno = unique(list_datasc_pro_SCType$test_neuron_auto@meta.data$coarse_cell_type))

# build the scatter plot
df_jaccard_score <- pmap(list(SCType_anno = df_crossing$SCType_anno,
                              cellranger_anno = df_crossing$cellranger_anno), function(SCType_anno,cellranger_anno){
                                
                                # calculate the jaccard score
                                a <- df_meta_full2 %>%
                                  filter(annotation_confident == SCType_anno) %>% pull(barcode)
                                b <- df_meta_full2 %>%
                                  filter(coarse_cell_type == cellranger_anno) %>% pull(barcode)
                                jaccard_score <- jaccard(a,b)
                                
                                # build a data.frame
                                df <- data.frame(SCType_anno = SCType_anno,
                                                 cellranger_anno = cellranger_anno,
                                                 jaccard_score = jaccard_score)
                                return(df)
                              }) %>%
  bind_rows()

head(df_jaccard_score)

# shape it as a matrix
mat_jaccard_score <- df_jaccard_score %>%
  pivot_wider(names_from = cellranger_anno,values_from = jaccard_score) %>%
  column_to_rownames("SCType_anno")

mat_jaccard_score

# plot the matrix
ht_02 <- Heatmap(mat_jaccard_score,
                 name = "Jaccard score",
                 # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                 col = viridis::viridis(option = "turbo",n = 20),
                 row_names_side = "right",
                 row_names_gp = gpar(fontsize = 8),
                 column_names_side = "bottom",
                 column_names_gp = gpar(fontsize = 8),
                 row_dend_reorder = FALSE,
                 column_dend_reorder = FALSE,
                 row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 show_column_names = T,
                 show_row_names = T)

ht_02
