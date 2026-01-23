# AIM ---------------------------------------------------------------------
# before integration compare the two object with or without introns

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
library(DESeq2)
library(edgeR)
library(RNAseqQC)
library(vegan)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
# load the LUT
# LUT <- read_csv("../../data/LUT_samples.csv")

# read in the list of objects. use the filtered dataset for the singlets only
data.list <- readRDS("../out/test_introns/object/list_datasc_fix_filter_norm_doubletSinglet_SoupX_01000_06000_15_V5.rds")

# EDA ---------------------------------------------------------------------

# compare metrics =========================================================
# Compare the number of features per cell. use the cell barcode as link
df_meta_full <- inner_join(
  data.list$sample_untreated_WOintron@meta.data %>%
    rownames_to_column("barcode") %>%
    dplyr::select(barcode,nFeature_RNA,percent.mt,percent.ribo,percent.globin),
  data.list$sample_untreated_Wintron@meta.data %>%
    rownames_to_column("barcode") %>%
    dplyr::select(barcode,nFeature_RNA,percent.mt,percent.ribo,percent.globin),suffix = c("_WOintron","_Wintron"),by = "barcode") %>%
  pivot_longer(names_to = "ID",values_to = "value",-barcode) %>%
  mutate(feature_name = str_extract(ID,"nFeature_RNA|percent.mt|percent.ribo|percent.globin"),
         test = str_extract(ID,"WOintron|Wintron"))

# plot violins per test
df_meta_full %>%
  ggplot(aes(x = test,y = value)) +
  geom_violin() +
  geom_boxplot(width = 0.1,outlier.shape = NA) +
  facet_wrap(~feature_name,scales = "free_y",nrow = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Features per cell",
       x = "Feature",
       y = "Value") +
  scale_y_sqrt()
ggsave("../out/test_introns/plot/BoxplotFeatures_before_integration_V5.pdf", width = 12, height = 3)

# plot the scatters
df_meta_full %>%
  select(-ID) %>%
  pivot_wider(names_from = test,values_from = value) %>%
  ggplot(aes(x = WOintron,y = Wintron)) +
  geom_point(size=0.1,alpha = 0.1) +
  facet_wrap(~feature_name,scales = "free") +
  theme_minimal() +
  theme(strip.background = element_blank()) +
  scale_y_sqrt()+
  scale_x_sqrt()+
  geom_abline(slope = 1,intercept = 0,linetype = "dashed",col="red")
ggsave("../out/test_introns/plot/ScatterFeatures_before_integration_V5.pdf", width = 7, height = 7)

# plot UMAPS =============================================================
# plot the two UMAPs before integration
p1 <- DimPlot(object = data.list$sample_untreated_WOintron, group.by = "seurat_clusters", label = T,raster = T) + NoLegend()+ ggtitle("Without introns")
p2 <- DimPlot(object = data.list$sample_untreated_Wintron, group.by = "seurat_clusters", label = T,raster = T) + NoLegend()+ ggtitle("With introns")

p1 + p2
ggsave("../out/test_introns/plot/UMAP_before_integration_V5.png", width = 10, height = 5)

# compare the clusters ===================================================
# build the matrix for the barcode assignament
df_meta_full2 <- inner_join(
  data.list$sample_untreated_WOintron@meta.data %>%
    rownames_to_column("barcode") %>%
    dplyr::select(barcode,seurat_clusters) %>%
    mutate(seurat_clusters = paste0("WOintron_",seurat_clusters)),
  data.list$sample_untreated_Wintron@meta.data %>%
    rownames_to_column("barcode") %>%
    dplyr::select(barcode,seurat_clusters)%>%
    mutate(seurat_clusters = paste0("Wintron_",seurat_clusters)),
  suffix = c("_WOintron","_Wintron"),by = "barcode")

# compare counts per cluster #############################################
# try to summarise the counts per cluster and shape it as a matrix
mat_counts <- df_meta_full2 %>%
  group_by(seurat_clusters_WOintron,seurat_clusters_Wintron) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = seurat_clusters_Wintron,values_from = n,values_fill = 0) %>%
  column_to_rownames("seurat_clusters_WOintron")

# t(scale(t(mat_counts))) %>% rowSums()
# t(scale(t(mat_counts))) %>% rowSds()

# for graphical purposes, scale the matrix by rows
ht_01 <- Heatmap(t(scale(t(mat_counts))),
              name = "scaled cell counts",
              # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
              col = viridis::viridis(option = "turbo",n = 50),
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

# Jaccard similarity #####################################################
# Jaccard score is defined as:
# Jaccard Similarity = (number of observations in both sets) / (number in either set)
# define the jaccard score function
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}
# # test
# a <- c('potato', 'tomotto', 'chips', 'baloon')
# b <- c('car', 'chips', 'bird', 'salt')
# jaccard(a, b)

# build all the comparison
# build the dataset for the correlatino plot
df_crossing <- crossing(seurat_clusters_WOintron = unique(df_meta_full2$seurat_clusters_WOintron),
                        seurat_clusters_Wintron = unique(df_meta_full2$seurat_clusters_Wintron))

# build the scatter plot
df_jaccard_score <- pmap(list(WOintron = df_crossing$seurat_clusters_WOintron,
                              Wintron = df_crossing$seurat_clusters_Wintron), function(WOintron,Wintron){
                                
                                # calculate the jaccard score
                                a <- df_meta_full2 %>%
                                  filter(seurat_clusters_WOintron == WOintron) %>% pull(barcode)
                                b <- df_meta_full2 %>%
                                  filter(seurat_clusters_Wintron == Wintron) %>% pull(barcode)
                                jaccard_score <- jaccard(a,b)
                                
                                # build a data.frame
                                df <- data.frame(WOintron = WOintron,
                                                 Wintron = Wintron,
                                                 jaccard_score = jaccard_score)
                                return(df)
                                }) %>%
  bind_rows()

# shape it as a matrix
mat_jaccard_score <- df_jaccard_score %>%
  pivot_wider(names_from = Wintron,values_from = jaccard_score) %>%
  column_to_rownames("WOintron")

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
pdf("../out/test_introns/plot/Heatmap_Jaccard_score_before_integration_V5.pdf", width = 7, height = 6)
draw(ht_02)
dev.off()

# correlation matrix =====================================================
# Build pseudobulk per clusters and attempt a correlation analysis cluster-wise
# aggregate the expression per sample per treatment, donor and cell type
mat_WOintron <- AggregateExpression(object = data.list$sample_untreated_WOintron,
                                    group.by = c("seurat_clusters"),
                                    assays = 'RNA',
                                    slot = "counts",
                                    return.seurat = FALSE)

mat_Wintron <- AggregateExpression(object = data.list$sample_untreated_Wintron,
                                    group.by = c("seurat_clusters"),
                                    assays = 'RNA',
                                    slot = "counts",
                                    return.seurat = FALSE)

# 1. Get counts matrix
counts_WOintron <- mat_WOintron$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene,names_to = "sample",values_to = "counts") %>%
  # rename the samples
  mutate(sample_new = paste0("WOintron_",sample)) %>%
  select(-sample) %>%
  pivot_wider(names_from = sample_new,values_from = counts) %>%
  column_to_rownames("gene")

counts_Wintron <- mat_Wintron$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene,names_to = "sample",values_to = "counts") %>%
  # rename the samples
  mutate(sample_new = paste0("Wintron_",sample)) %>%
  select(-sample) %>%
  pivot_wider(names_from = sample_new,values_from = counts) %>%
  column_to_rownames("gene")

# 2. generate sample level metadata
LUT_WOintron <- data.frame(sample = colnames(counts_WOintron))
LUT_Wintron <- data.frame(sample = colnames(counts_Wintron))

# 3. build the deseq2 object
# build the design matrix
design_WOintron <- model.matrix(~ LUT_WOintron$sample)
colnames(design_WOintron)[1] <- c("intercept")

design_Wintron <- model.matrix(~ LUT_Wintron$sample)
colnames(design_Wintron)[1] <- c("intercept")

# build the object
dds_WOintron <- DESeqDataSetFromMatrix(countData = counts_WOintron,
                                       colData = LUT_WOintron,
                                       design = design_WOintron)

dds_Wintron <- DESeqDataSetFromMatrix(countData = counts_Wintron,
                                      colData = LUT_Wintron,
                                      design = design_Wintron)

# 4. remove lowly expressed genes
keep_WOintron <- edgeR::filterByExpr(counts(dds_WOintron), group = LUT_WOintron$sample)
dds_WOintron_filter <- dds_WOintron[keep_WOintron,]

keep_Wintron <- edgeR::filterByExpr(counts(dds_Wintron), group = LUT_Wintron$sample)
dds_Wintron_filter <- dds_Wintron[keep_Wintron,]

# plot the raw number of reads per sample
colSums(counts(dds_WOintron_filter)) %>%
  data.frame(tot_reads = .) %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x=sample,y=tot_reads)) + geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))

colSums(counts(dds_Wintron_filter)) %>%
  data.frame(tot_reads = .) %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x=sample,y=tot_reads)) + geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))

# 5. scale the data
vds_WOintron_filter <- vst(dds_WOintron_filter, blind = T)
vds_Wintron_filter <- vst(dds_Wintron_filter, blind = T)

# on the scaled data run the correlation using the common genes
# define a common set of genes to make the correlation matrix
common_genes <- intersect(rownames(assay(vds_WOintron_filter)),rownames(assay(vds_Wintron_filter)))

# generate the two matrices as ordered by the same set of features
mat_WOintron_filter <- assay(vds_WOintron_filter)[common_genes,]
mat_Wintron_filter <- assay(vds_Wintron_filter)[common_genes,]

# make sure both the dataset have the same rownames
sum(!rownames(mat_WOintron_filter) == rownames(mat_Wintron_filter))

# build the correlation matrix per sample
cor_mat <- cor(mat_WOintron_filter,mat_Wintron_filter)

ht_03 <- Heatmap(cor_mat,
              name = "Correlation \nPearson",
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

# Attempted annotation based on marker genes =============================
shortlist_features_list2 <- list(
  IMMUNE = c("AIF1","TYROBP","HLA-DRA","TREM2","CX3CR1","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("CSPG4","OLIG1","OLIG2", "PDGFRA", "SOX6", "PLP1","SOX10","BCAS1","MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "VIM","SLC1A2","S100B"),
  # Neu = c("SYT1")
  NEURONS = c("GAD2", "TLE4", "CUX2","SYP", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF")
)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_short_WOintron <- DotPlot(data.list$sample_untreated_WOintron, features = shortlist_features_list2, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()+ ggtitle("Without introns")
test_short_Wintron <- DotPlot(data.list$sample_untreated_Wintron, features = shortlist_features_list2, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()+ ggtitle("With introns")

test_short_WOintron/test_short_Wintron
ggsave("../out/test_introns/plot/Dotplot_markers_before_integration_V5.pdf",width = 13,height = 10)

# based on the expression of the markers gene try to annotate the clusters
# IMPORTANT: make sure not to lose the rowname in the metadata
data.list$sample_untreated_WOintron@meta.data <-data.list$sample_untreated_WOintron@meta.data %>%
  mutate(cell_id = case_when(seurat_clusters %in% c(6,2,11,13)~"ASTRO",
                             seurat_clusters %in% c(9,5,10,14,3)~"NEU",
                             seurat_clusters %in% c(0)~"OPC",
                             seurat_clusters %in% c(8)~"PROG",
                             seurat_clusters %in% c(4)~"OLIGO",
                             seurat_clusters %in% c(12)~"MG",
                             seurat_clusters %in% c(1,7,16)~"GLIA_IMM"))

data.list$sample_untreated_Wintron@meta.data <- data.list$sample_untreated_Wintron@meta.data %>%
  mutate(cell_id = case_when(seurat_clusters %in% c(6,3,11,15,19)~"ASTRO",
                             seurat_clusters %in% c(2,5,9,10,14,17,18)~"NEU",
                             seurat_clusters %in% c(0)~"OPC",
                             seurat_clusters %in% c(8)~"PROG",
                             seurat_clusters %in% c(4)~"OLIGO",
                             seurat_clusters %in% c(12)~"MG",
                             seurat_clusters %in% c(1,7,16,20)~"GLIA_IMM"))

# plot the two UMAPs before integration
p12 <- DimPlot(object = data.list$sample_untreated_WOintron, group.by = "cell_id", label = T,raster = T) + NoLegend()+ ggtitle("Without introns")
p22 <- DimPlot(object = data.list$sample_untreated_Wintron, group.by = "cell_id", label = T,raster = T) + NoLegend()+ ggtitle("With introns")

p12 + p22
ggsave("../out/test_introns/plot/UMAP_before_integration_annotatation_V5.png", width = 10, height = 5)

# Jaccard similarity #####################################################
# build the matrix for the barcode assignament
df_meta_full2_annotation <- inner_join(
  data.list$sample_untreated_WOintron@meta.data %>%
    rownames_to_column("barcode") %>%
    dplyr::select(barcode,cell_id) %>%
    mutate(seurat_clusters = paste0("WOintron_",cell_id)),
  data.list$sample_untreated_Wintron@meta.data %>%
    rownames_to_column("barcode") %>%
    dplyr::select(barcode,cell_id)%>%
    mutate(seurat_clusters = paste0("Wintron_",cell_id)),
  suffix = c("_WOintron","_Wintron"),by = "barcode")

df_meta_full2_annotation %>%
  filter(cell_id_WOintron != cell_id_Wintron)

# Jaccard score is defined as:
# Jaccard Similarity = (number of observations in both sets) / (number in either set)
# define the jaccard score function
# jaccard <- function(a, b) {
#   intersection <- length(intersect(a, b))
#   union <- length(a) + length(b) - intersection
#   return (intersection/union)
# }
# # test
# a <- c('potato', 'tomotto', 'chips', 'baloon')
# b <- c('car', 'chips', 'bird', 'salt')
# jaccard(a, b)

# build all the comparison
# build the dataset for the correlatino plot
df_crossing_annotation <- crossing(seurat_clusters_WOintron = unique(df_meta_full2_annotation$seurat_clusters_WOintron),
                                   seurat_clusters_Wintron = unique(df_meta_full2_annotation$seurat_clusters_Wintron))

# build the scatter plot
df_jaccard_score_annotation <- pmap(list(WOintron = df_crossing_annotation$seurat_clusters_WOintron,
                                         Wintron = df_crossing_annotation$seurat_clusters_Wintron), function(WOintron,Wintron){
                                
                                # calculate the jaccard score
                                a <- df_meta_full2_annotation %>%
                                  filter(seurat_clusters_WOintron == WOintron) %>% pull(barcode)
                                b <- df_meta_full2_annotation %>%
                                  filter(seurat_clusters_Wintron == Wintron) %>% pull(barcode)
                                jaccard_score <- jaccard(a,b)
                                
                                # build a data.frame
                                df <- data.frame(WOintron = WOintron,
                                                 Wintron = Wintron,
                                                 jaccard_score = jaccard_score)
                                return(df)
                              }) %>%
  bind_rows()

# shape it as a matrix
mat_jaccard_score_annotation <- df_jaccard_score_annotation %>%
  pivot_wider(names_from = Wintron,values_from = jaccard_score) %>%
  column_to_rownames("WOintron")

ht_02_annotation <- Heatmap(mat_jaccard_score_annotation,
                 name = "Jaccard score",
                 # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                 col = viridis::viridis(option = "turbo",n = 20,direction = -1),
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
pdf("../out/test_introns/plot/Heatmap_Jaccard_score_before_integration_V5_annotationInverted.pdf", width = 5, height = 4)
draw(ht_02_annotation)
dev.off()

png("../out/test_introns/plot/Heatmap_Jaccard_score_before_integration_V5_annotationInverted.png",width = 500,height = 400)
draw(ht_02_annotation)
dev.off()

# -------------------------------------------------------------------------
# dotplot after marker anntoation
test_short_WOintron2 <- DotPlot(data.list$sample_untreated_WOintron, features = shortlist_features_list2, dot.scale = 8,cluster.idents = T,group.by = "cell_id") +
  RotatedAxis()+ ggtitle("Without introns")
test_short_Wintron2 <- DotPlot(data.list$sample_untreated_Wintron, features = shortlist_features_list2, dot.scale = 8,cluster.idents = T,group.by = "cell_id") +
  RotatedAxis()+ ggtitle("With introns")

test_short_WOintron2/test_short_Wintron2
ggsave("../out/test_introns/plot/Dotplot_markers_before_integration_annotation_V5.pdf",width = 13,height = 10)

# build the matrix for the barcode assignament
df_meta_full3 <- inner_join(
  data.list$sample_untreated_WOintron@meta.data %>%
    rownames_to_column("barcode") %>%
    dplyr::select(barcode,cell_id) %>%
    mutate(cell_id = paste0("WOintron_",cell_id)),
  data.list$sample_untreated_Wintron@meta.data %>%
    rownames_to_column("barcode") %>%
    dplyr::select(barcode,cell_id) %>%
    mutate(cell_id = paste0("Wintron_",cell_id)),
  suffix = c("_WOintron","_Wintron"),by = "barcode")

# compare counts per cluster #############################################
# try to summarise the counts per cluster and shape it as a matrix
mat_counts3 <- df_meta_full3 %>%
  group_by(cell_id_WOintron,cell_id_Wintron) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = cell_id_Wintron,values_from = n,values_fill = 0) %>%
  column_to_rownames("cell_id_WOintron")

# t(scale(t(mat_counts))) %>% rowSums()
# t(scale(t(mat_counts))) %>% rowSds()

# for graphical purposes, scale the matrix by rows
ht_04 <- Heatmap(t(scale(t(mat_counts3))),
                 name = "scaled cell counts",
                 # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                 col = viridis::viridis(option = "turbo",n = 50),
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
pdf("../out/test_introns/plot/Heatmap_scaledCount_score_before_integration_annotation_V5.pdf", width = 5, height = 4)
draw(ht_04)
dev.off()

# Jaccard similarity #####################################################
# build all the comparison
# build the dataset for the correlatino plot
df_crossing2 <- crossing(cell_id_WOintron = unique(df_meta_full3$cell_id_WOintron),
                         cell_id_Wintron = unique(df_meta_full3$cell_id_Wintron))

# build the scatter plot
df_jaccard_score2 <- pmap(list(WOintron = df_crossing2$cell_id_WOintron,
                               Wintron = df_crossing2$cell_id_Wintron), function(WOintron,Wintron){
                                
                                # calculate the jaccard score
                                a <- df_meta_full3 %>%
                                  filter(cell_id_WOintron == WOintron) %>% pull(barcode)
                                b <- df_meta_full3 %>%
                                  filter(cell_id_Wintron == Wintron) %>% pull(barcode)
                                jaccard_score <- jaccard(a,b)
                                
                                # build a data.frame
                                df <- data.frame(WOintron = WOintron,
                                                 Wintron = Wintron,
                                                 jaccard_score = jaccard_score)
                                return(df)
                              }) %>%
  bind_rows()

# shape it as a matrix
mat_jaccard_score2 <- df_jaccard_score2 %>%
  pivot_wider(names_from = Wintron,values_from = jaccard_score) %>%
  column_to_rownames("WOintron")

ht_05 <- Heatmap(mat_jaccard_score2,
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
pdf("../out/test_introns/plot/Heatmap_Jaccard_score_before_integration_annotation_V5.pdf", width = 7, height = 6)
draw(ht_05)
dev.off()

# Identify conserved cell type markers ------------------------------------
# data
DefaultAssay(data.list$sample_untreated_WOintron) <- "RNA"
DefaultAssay(data.list$sample_untreated_Wintron) <- "RNA"

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Idents(data.list$sample_untreated_WOintron) <- "cell_id"
sobj_total_h.markers_WOintron <- RunPrestoAll(data.list$sample_untreated_WOintron, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Idents(data.list$sample_untreated_Wintron) <- "cell_id"
sobj_total_h.markers_Wintron <- RunPrestoAll(data.list$sample_untreated_Wintron, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# compare the tables per annotation
# number of marker genes per cluster
list(sobj_total_h.markers_WOintron %>%
       select(cluster,gene) %>%
       mutate(test = "WOintron"),
     sobj_total_h.markers_Wintron %>%
       select(cluster,gene) %>%
       mutate(test = "Wintron")) %>%
  bind_rows() %>%
  group_by(cluster,test) %>% summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(cluster = fct_reorder(cluster,-tot)) %>%
  ggplot(aes(x=cluster,y=n,fill = test)) + geom_col(position = "dodge")+theme_bw()+
  labs(title = "number of marker genes per cell id (RunPrestoAll)",
       x = "cell_id",
       y = "number of marker gene")
ggsave("../out/test_introns/plot/Barplot_marker_genes_before_integration_annotation_V5.pdf",width = 10,height = 5)

# comapre the identity of the marker genes
cell_id <- unique(sobj_total_h.markers_WOintron$cluster)

df_gene_comparison <- lapply(cell_id, function(x){
  gene_WOintron <- sobj_total_h.markers_WOintron %>%
    filter(cluster %in% x) %>%
    pull(gene)
  
  gene_Wintron <- sobj_total_h.markers_Wintron %>%
    filter(cluster %in% x) %>%
    pull(gene)
  
  df <- data.frame(common = length(intersect(gene_WOintron,gene_Wintron)),
                   WOintron_only = length(setdiff(gene_WOintron,gene_Wintron)),
                   Wintron_only = length(setdiff(gene_Wintron,gene_WOintron)),
                   cell_id = x)
    
}) %>%
  bind_rows()

df_gene_comparison %>%
  pivot_longer(names_to = "intersection",values_to = "value",-cell_id) %>%
  group_by(cell_id) %>%
  mutate(tot = sum(value)) %>%
  ungroup() %>%
  mutate(cell_id = fct_reorder(cell_id,-tot)) %>%
  ggplot(aes(x=cell_id,y=value,fill = intersection)) + geom_col(position = "dodge")+
  theme_bw()+
  labs(title = "number of marker genes per cell id (RunPrestoAll)",
       x = "cell_id",
       y = "number of marker gene")+
  scale_fill_manual(values = c("magenta","blue","red"))
ggsave("../out/test_introns/plot/BarplotCommon_marker_genes_before_integration_annotation_V5.pdf",width = 10,height = 5)

# check the reson behind the WOintron only
Idents(data.list$sample_untreated_WOintron) <- "cell_id"
sobj_total_h.markers_WOintron2 <- RunPrestoAll(data.list$sample_untreated_WOintron, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0)
Idents(data.list$sample_untreated_Wintron) <- "cell_id"
sobj_total_h.markers_Wintron2 <- RunPrestoAll(data.list$sample_untreated_Wintron, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0)


GOI_WOintron <- lapply(cell_id, function(x){
  gene_WOintron <- sobj_total_h.markers_WOintron %>%
    filter(cluster %in% x) %>%
    pull(gene)
  
  gene_Wintron <- sobj_total_h.markers_Wintron %>%
    filter(cluster %in% x) %>%
    pull(gene)
  
  df <- data.frame(gene = setdiff(gene_WOintron,gene_Wintron),
                   cluster = x) %>%
    left_join(sobj_total_h.markers_Wintron2 %>% select(gene,cluster,avg_log2FC,p_val_adj),by = c("gene","cluster")) %>%
    left_join(sobj_total_h.markers_WOintron2 %>% select(gene,cluster,avg_log2FC,p_val_adj),by = c("gene","cluster"),suffix=c(".Wintron",".WOintron"))
  
}) %>%
  bind_rows()

GOI_WOintron %>%
  ggplot(aes(x=avg_log2FC.Wintron,y=avg_log2FC.WOintron)) + geom_point() + facet_wrap(~cluster)

GOI_WOintron %>%
  ggplot(aes(x=p_val_adj.Wintron,y=p_val_adj.WOintron)) + geom_point() + facet_wrap(~cluster)
