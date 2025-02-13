# AIM ---------------------------------------------------------------------
# sample processing of Azimuth processing for reference labelling
# this comes after running the sample_Azimuth_01_preprocessing.R script
# this part run the actual label tranfer

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
options(Seurat.object.assay.version = "v3")

# run azimuth -------------------------------------------------------------
query <- LoadFileInput("../out/test_introns/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15_V3_diet.rds")
head(query@meta.data)

# The RunAzimuth function can take a Seurat object as input
# pass the folder with the reference data
# Name of reference to map to or a path to a directory containing ref.Rds and idx.annoy
query <- RunAzimuth(query, reference = "renv/library/R-4.3/x86_64-pc-linux-gnu/humancortexref.SeuratData/azimuth/")

# can also be passed
# query <- RunAzimuth(query, reference = "humancortexref")
saveRDS(query,"../out/object/Azimuth_humancortexref_harmonySkipIntegration_AllSoupX_01000_06000_15.rds")

# the metadata have all the info of interest
head(query@meta.data)

DimPlot(query, group.by = "predicted.subclass", label = TRUE, label.size = 3,raster = T) + NoLegend()
DimPlot(query, group.by = "predicted.subclass", label = TRUE, label.size = 3,reduction = "ref.umap",raster = T) + NoLegend()
DimPlot(query, group.by = "predicted.subclass", label = TRUE, label.size = 3,split.by = "orig.ident",raster = T) + NoLegend()

# save the UMAP coordinates of the new reference
query@reductions$ref.umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column("barcode") %>% 
  write_tsv("../out/table/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef.tsv")

# save the new annotation from azimuth
query@meta.data %>% 
  data.frame() %>% 
  rownames_to_column("barcode") %>% 
  write_tsv("../out/table/harmonySkipIntegration_AllSoupX_01000_06000_15_meta_AzimuthHumanCortexRef.tsv")

# costume plots threshold 0.75 --------------------------------------------
# costume visualization
UMAP_query <- read_tsv("../out/table/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef.tsv")
meta_query <- read_tsv("../out/table/harmonySkipIntegration_AllSoupX_01000_06000_15_meta_AzimuthHumanCortexRef.tsv") %>%
  mutate(robust_score_subclass = case_when(predicted.subclass.score>0.75&mapping.score>0.75~predicted.subclass,
                                           T~"uncertain"),
         robust_score_class = case_when(predicted.class.score>0.75&mapping.score>0.75~predicted.class,
                                        T~"uncertain"),
         robust_score_cluster = case_when(predicted.cluster.score>0.75&mapping.score>0.75~predicted.cluster,
                                          T~"uncertain"))

# define the levels of the cluster variable
level_subclass <- meta_query %>%
  group_by(predicted.subclass) %>%
  summarise(med = median(predicted.subclass.score)) %>%
  mutate(predicted.subclass = fct_reorder(predicted.subclass, med,.desc = T)) %>%
  pull(predicted.subclass) %>%
  levels()

# for each assignment what is the distribution of the scores
meta_query %>%
  mutate(predicted.subclass = factor(predicted.subclass, levels = level_subclass)) %>%
  ggplot(aes(x=predicted.subclass,y=predicted.subclass.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
  geom_hline(yintercept = 0.75,col="red")
ggsave("../out/plot/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef_subclass_score_075.pdf",height = 4,width = 4)

meta_query %>%
  mutate(predicted.subclass = factor(predicted.subclass, levels = level_subclass)) %>%
  ggplot(aes(y=predicted.subclass,x=predicted.subclass.score))+
  ggridges::geom_density_ridges()+
  theme_bw() + 
  geom_vline(xintercept = 0.75,col="red",linetype="dashed")
ggsave("../out/plot/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef_subclass_ridges_score_075.pdf",height = 5,width = 4)

# define the levels of the cluster variable
level_class <- meta_query %>%
  group_by(predicted.class) %>%
  summarise(med = median(predicted.class.score)) %>%
  mutate(predicted.class = fct_reorder(predicted.class, med,.desc = T)) %>%
  pull(predicted.class) %>%
  levels()

# plot the score varaible for the cluster
meta_query %>%
  mutate(predicted.class = factor(predicted.class, levels = level_class)) %>%
  ggplot(aes(x=predicted.class,y=predicted.class.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
  geom_hline(yintercept = 0.75,col="red")
ggsave("../out/plot/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef_class_score_075.pdf",height = 4,width = 2)

# define the levels of the cluster variable
level_cluster <- meta_query %>%
  group_by(predicted.cluster) %>%
  summarise(med = median(predicted.cluster.score)) %>%
  mutate(predicted.cluster = fct_reorder(predicted.cluster, med,.desc = T)) %>%
  pull(predicted.cluster) %>%
  levels()

# plot the score varaible for the cluster
meta_query %>%
  mutate(predicted.cluster = factor(predicted.cluster, levels = level_cluster)) %>%
  ggplot(aes(x=predicted.cluster,y=predicted.cluster.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45),
        plot.margin = margin(0.5, 0.5, 0.5, 2, "cm")) + 
  geom_hline(yintercept = 0.75,col="red")
ggsave("../out/plot/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef_cluster_score_075.pdf",height = 4,width = 6)

# meta_query %>%
#   ggplot(aes(x=predicted.subclass,y=predicted.subclass.score))+
#   geom_violin(scale = "width") +
#   theme_bw() + 
#   theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
#   geom_hline(yintercept = 0.7,col="red")

# identifyt he most likely assignment for each seurat cluster
# first using all the subcluster annotation, not filtering for threshold of scores
prop_table_subclass <- meta_query %>%
  group_by(seurat_clusters,predicted.subclass) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.subclass,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.subclass") %>% 
  as.matrix()

pdf("../out/plot/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef_subclass_Heatmap_075.pdf",height = 4,width = 4)
Heatmap(prop_table_subclass,
        name = "prop", 
        column_title = "subclass score all",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

meta_query %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters,predicted.subclass) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

# first using all the subcluster annotation, filtering for threshold of scores
prop_table_subclass_filter <- meta_query %>%
  # filter(robust_score_subclass != "uncertain") %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters,predicted.subclass) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.subclass,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.subclass") %>%
  as.matrix()

pdf("../out/plot/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef_subclass_HeatmapFilter_075.pdf",height = 4,width = 5)
Heatmap(prop_table_subclass_filter,
        name = "prop", 
        column_title = "subclass score high confidence",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

meta_query %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters,predicted.subclass) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.subclass,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.subclass")

# add the meta to the coordinates
data_query <- left_join(UMAP_query,meta_query,"barcode")
# add the groupign suggested by Martina
# data_query$group_martina <- case_when(data_query$orig.ident %in% c("pool_MOCK_CTRL")~"MOCK",
#                                       data_query$orig.ident %in% c("CTRL8_WT_CSF","CTRL8_WT_CTRL")~"CTRL",
#                                       T~"SOX10")

# divide the dataset into uncertain and not
data_query_unc <- data_query %>%
  filter(robust_score_subclass == "uncertain")
#
data_query_defined <- data_query %>%
  filter(robust_score_subclass != "uncertain")

# average the position of the clusters
data_query_avg <- data_query_defined %>% group_by(robust_score_subclass) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)
data_query_avg2 <- data_query_defined %>% group_by(seurat_clusters) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score_subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg,aes(x = UMAP_1,y = UMAP_2,label = robust_score_subclass),col="black")+theme_bw()
ggsave("../out/plot/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef_subclass_UMAP_075.pdf",height = 5,width = 3)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score_subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg,aes(x = UMAP_1,y = UMAP_2,label = robust_score_subclass),col="black")+theme_bw()+
  facet_wrap(~ID,ncol=4)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../out/plot/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef_subclass_split_UMAP_075.pdf",height = 15,width = 20)

# data_query %>%
#   mutate(seurat_clusters=factor(seurat_clusters)) %>%
#   ggplot(label= TRUE)+
#   # geom_point(data = data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
#   geom_point(aes(x = UMAP_1,y = UMAP_2, col = seurat_clusters),size=0.3,alpha=0.8) +
#   guides(colour = guide_legend(override.aes = list(size=5))) +
#   # labs(color= "Clusters") +
#   ggrepel::geom_text_repel(data = data_query_avg2,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+theme_bw()
# ggsave("out/image/data.combined_annotated_norm_fix_regressCC_resolution_DoubletSinglet_UMAP.pdf",width = 8,height = 6)

# costum plotting ---------------------------------------------------------
# costume visualization
UMAP_ref <- read_tsv("../out/test_introns/table/azimuth_humancortex_CoordUMAP.tsv")
UMAP_ref2 <- read_tsv("../out/table/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef.tsv")
meta_ref <- read_tsv("../out/test_introns/table/azimuth_humancortex_Metadata.tsv")
meta_query <- read_tsv("../out/table/harmonySkipIntegration_AllSoupX_01000_06000_15_meta_AzimuthHumanCortexRef.tsv") %>%
  mutate(robust_score_subclass = case_when(predicted.subclass.score>0.75&mapping.score>0.75~predicted.subclass,
                                           T~"uncertain"),
         robust_score_class = case_when(predicted.class.score>0.75&mapping.score>0.75~predicted.class,
                                        T~"uncertain"),
         robust_score_cluster = case_when(predicted.cluster.score>0.75&mapping.score>0.75~predicted.cluster,
                                          T~"uncertain"))

# add the meta to the coordinates
data_ref <- left_join(UMAP_ref,meta_ref,"barcode")
data_ref2 <- left_join(UMAP_ref2,meta_query,"barcode") %>%
  mutate(seurat_clusters=factor(seurat_clusters))

# average the position of the clusters
data_ref_avg <- data_ref %>% group_by(subclass) %>% select(refumap_1,refumap_2) %>% summarize_all(mean)

# average the position of the clusters
data_ref_avg2 <- data_ref2 %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# plot the reference
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2,col=subclass),size=0.3,alpha=1) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg,aes(x = refumap_1,y = refumap_2,label = subclass),col="black",max.overlaps = 12)+theme_cowplot()
ggsave("../out/plot/human_motorcortexv1.0.0_subclass_alt.pdf",width = 8,height = 6)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = seurat_clusters),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black",max.overlaps = 12)+theme_bw()
ggsave("../out/plot/human_motorcortexv1.0.0_subclass_data.combined_annotated_norm_fix_resolution_UMAP_075.pdf",width = 9,height = 6)

ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = seurat_clusters),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black",max.overlaps = 15)+theme_bw()+facet_wrap(~ID,ncol=4)+theme(strip.background = element_blank())
ggsave("../out/plot/human_motorcortexv1.0.0_subclass_data.combined_annotated_norm_fix_resolution_UMAP_075_split.pdf",width = 20,height = 15)

# average the position of the clusters
data_ref_avg2subclass <- data_ref2 %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(predicted.subclass) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2subclass,aes(x = UMAP_1,y = UMAP_2,label = predicted.subclass),col="black")+theme_bw()
ggsave("../out/plot/human_motorcortexv1.0.0_subclass_data.combined_annotated_norm_fix_resolution_UMAP_075_subclass.pdf",width = 8,height = 6)

ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2subclass,aes(x = UMAP_1,y = UMAP_2,label = predicted.subclass),col="black")+theme_bw()+facet_wrap(~ID,ncol=4)+theme(strip.background = element_blank())
ggsave("../out/plot/human_motorcortexv1.0.0_subclass_data.combined_annotated_norm_fix_resolution_UMAP_075_subclass_split.pdf",width = 20,height = 15)

# inverse transfer, see the ref on the query ------------------------------
# plot the annotation over the original data
# read in the original dataset
# "../out/test_introns/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15_V3_diet.rds"
sobj <- readRDS("../out/test_introns/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15_V3.rds")

# load the meta with annotations
meta_query <- read_tsv("../out/table/harmonySkipIntegration_AllSoupX_01000_06000_15_meta_AzimuthHumanCortexRef.tsv") %>%
  mutate(robust_score_subclass = case_when(predicted.subclass.score>0.75&mapping.score>0.75~predicted.subclass,
                                           T~"uncertain"),
         robust_score_class = case_when(predicted.class.score>0.75&mapping.score>0.75~predicted.class,
                                        T~"uncertain"),
         robust_score_cluster = case_when(predicted.cluster.score>0.75&mapping.score>0.75~predicted.cluster,
                                          T~"uncertain"))

test_UMAP_query <- sobj@reductions$umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column("barcode")

test_data_query <- left_join(test_UMAP_query,meta_query,"barcode")
# add the groupign suggested by Martina
# data_query$group_martina <- case_when(data_query$orig.ident %in% c("pool_MOCK_CTRL")~"MOCK",
#                                       data_query$orig.ident %in% c("CTRL8_WT_CSF","CTRL8_WT_CTRL")~"CTRL",
#                                       T~"SOX10")

# divide the dataset into uncertain and not
test_data_query_unc <- test_data_query %>%
  filter(robust_score_subclass == "uncertain")

#
test_data_query_defined <- test_data_query %>%
  filter(robust_score_subclass != "uncertain")

# average the position of the clusters
test_data_query_avg <- test_data_query_defined %>% group_by(robust_score_subclass) %>% select(umap_1, umap_2) %>% summarize_all(mean)
test_data_query_avg2 <- test_data_query_defined %>% group_by(seurat_clusters) %>% select(umap_1, umap_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = test_data_query_unc,aes(x = umap_1,y = umap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = test_data_query_defined,aes(x = umap_1,y = umap_2, col = robust_score_subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = test_data_query_avg,aes(x = umap_1,y = umap_2,label = robust_score_subclass),col="black")+theme_bw()
ggsave("../out/plot/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthBS_subclass_UMAP_075.pdf",height = 3,width = 5)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = test_data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = test_data_query_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score_subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = test_data_query_avg,aes(x = UMAP_1,y = UMAP_2,label = robust_score_subclass),col="black")+theme_bw()+
  facet_wrap(~ID,ncol=4)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthBS_subclass_split_UMAP_075.pdf",height = 15,width = 20)
