# AIM ---------------------------------------------------------------------
# sample routine for the label transfer approach

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(SeuratData)
library(ggridges)
library(ComplexHeatmap)
library(SeuratWrappers)
library(cowplot)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v3")

# read in the data --------------------------------------------------------
# read in a sample object to use as reference
ref_WM <- readRDS("../data/all20_immune_model.rds")

# check the reference
DimPlot(ref_WM,label = T,group.by = "cellID")

# read in the query dataset this is the clean3 version of the dataset
query <- readRDS("../data/129_MG_subcluster_HarmonySample.rds")

# check the query
DimPlot(query,label = T,raster = T,group.by = "RNA_snn_res.0.9")

# if needed fix the meta to add the general cell annotation
query$seurat_clusters_fix1 <- paste0("clu_",query$RNA_snn_res.0.9)
query@meta.data


# run the label transfer: ref over query ----------------------------------

# use the RNA slot to avoid issue in the FindTransferAnchors call
DefaultAssay(ref_WM)<-"RNA"

dim(ref_WM)
ref_WM@meta.data

# in the vignetted the query isn't the RNA slot
# DefaultAssay(query) <- "integrated"
# use the RNA slot in general to make sure all the features are accounted
DefaultAssay(query) <- "RNA"
dim(query)

# find the anchors
test.anchors_0 <- FindTransferAnchors(reference = ref_WM,
                                      # k.filter = NA,
                                      query = query,
                                      dims = 1:30,
                                      reference.reduction = "pca")

# predict the labels from the reference onto the query
predictions_0 <- TransferData(anchorset = test.anchors_0,
                              refdata = ref_WM$cellID,
                              dims = 1:30)

# notice that the dimension of the prediction is the same as the query object (in terms of number of cells)
dim(predictions_0)
dim(query)

# wrangling ---------------------------------------------------------------
# add the coordinates to the metadata
query_transfer <- AddMetaData(query, metadata = query@reductions$umap@cell.embeddings)

# add the predictions ot the query dataset
query_transfer <- AddMetaData(query_transfer, metadata = predictions_0)

# add the robust score metric to the object. use the 0.75 cut off as referenced in Azimuth https://azimuth.hubmapconsortium.org/
query_transfer$robust_score <- query_transfer@meta.data %>%
  rownames_to_column(var = "barcodes") %>%
  # if above the threshold, confirm the prediction, otherwise call uncertain
  mutate(robust_score = case_when(prediction.score.max>0.75~predicted.id,
                                  T ~ "uncertain")) %>%
  pull(robust_score)

# save the object
saveRDS(query_transfer,"../out/object/129_query_transfer_RefOverQuery_res0.9.rds")

# EDA of the tranferring --------------------------------------------------

# plot the predicition score =============================================
# show the distribution of the max scores
predictions_0 %>%
  ggplot(aes(x=prediction.score.max))+geom_histogram()+theme_bw()
# ggsave("../../out/image/129_histo_prediction_label_transfer_refWMImmune_queryIMMLYMsubset.pdf",width = 5,height = 4)

# check the table of the robust score
table(query_transfer$robust_score)

# default plot 
# plot the raw prediction based on the reference dataset
DimPlot(query_transfer,label = T,group.by = "predicted.id")

# plot the robust prediction
DimPlot(query_transfer,label = T,group.by = "robust_score")

# custom plot
# average the position of the clusters
data_transfer_avg_raw <- query_transfer@meta.data %>% group_by(predicted.id) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# plot the raw prediction based on the reference dataset
ggplot(label= TRUE)+
  geom_point(data = query_transfer@meta.data,aes(x = UMAP_1,y = UMAP_2, col = predicted.id),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_transfer_avg_raw,aes(x = UMAP_1,y = UMAP_2,label = predicted.id),col="black")+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5)))
# ggsave("../../out/image/129_UMAP_label_transfer_refWMImmune_queryIMMLYMsubset_raw.pdf",width = 5,height = 4)

# plot the robust prediction based on the reference dataset
# split the data for convinience
data_transfer_unc <- query_transfer@meta.data %>% 
  filter(robust_score == "uncertain")
#
data_transfer_defined <- query_transfer@meta.data %>% 
  filter(robust_score != "uncertain")

data_transfer_avg_robust <- data_transfer_defined %>% group_by(predicted.id) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

ggplot(label= TRUE)+
  # plot the uncertain only in the background as gray
  geom_point(data = data_transfer_unc,aes(x = UMAP_1,y = UMAP_2),size=0.1,alpha=0.1,col="gray") +
  geom_point(data = data_transfer_defined,aes(x = UMAP_1,y = UMAP_2, col = predicted.id),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_transfer_avg_robust,aes(x = UMAP_1,y = UMAP_2,label = predicted.id),col="black")+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5)))
# ggsave("../../out/image/129_UMAP_label_transfer_refWMImmune_queryIMMLYMsubset_robust.pdf",width = 5,height = 4)

# plot the predicted score per per cell type =============================
# since I am placing my cells on a reference it makes sense to me to show the score of how well each of my cells i mappe on the reference
id_factor <- query_transfer@meta.data %>% 
  group_by(seurat_clusters_fix1) %>% 
  summarise(med = median(prediction.score.max)) %>% 
  arrange(med) %>% 
  pull(seurat_clusters_fix1)

query_transfer@meta.data %>%
  mutate(seurat_clusters_fix1= factor(seurat_clusters_fix1,levels = id_factor)) %>% 
  ggplot(aes(x=prediction.score.max,y=seurat_clusters_fix1)) +
  geom_density_ridges(scale = 0.9) +
  theme_bw() +
  geom_vline(xintercept = 0.70,linetype="dashed",col="red",alpha=0.5)
# ggsave("../out/plot/129_histo_label_transfer_refWMImmune_queryIMMLYMsubset_res0.9.pdf",width = 3,height = 2)

# split by condition
query_transfer@meta.data %>%
  mutate(seurat_clusters_fix1= factor(seurat_clusters_fix1,levels = id_factor)) %>% 
  ggplot(aes(x=prediction.score.max,y=seurat_clusters_fix1)) +
  geom_density_ridges(scale = 0.9) +
  theme_bw() +
  geom_vline(xintercept = 0.70,linetype="dashed",col="red",alpha=0.5)+facet_wrap(~pathology_class)+theme(strip.background = element_blank())
# ggsave("../out/plot/129_histo_label_transfer_refWMImmune_queryIMMLYMsubset_split_res0.9.pdf",width = 9,height = 6)

# split the proportions per pathology_classment
prop_table_tot_test <- query_transfer@meta.data %>%
  # filter(robust_score_subclass != "uncertain") %>%
  # filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters_fix1,pathology_class,predicted.id,robust_score) %>%
  summarise(n = n()) %>%
  ungroup() %>% 
  group_by(seurat_clusters_fix1,pathology_class) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>% 
  ungroup() %>% 
  mutate(robust = case_when(robust_score!="uncertain"~"robust",
                            T~robust_score))

prop_table_tot_test %>% 
  filter(robust == "robust")

prop_table_tot_test %>% 
  ggplot(aes(x=pathology_class,y=prop,col=seurat_clusters_fix1))+
  geom_point(position = position_jitter(width = 0.2))+
  facet_grid(robust~predicted.id)+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45),strip.background = element_blank())
# ggsave("../out/plot/label_transfer_refWMImmune_queryBrainsphere_prop_table_tot_test.pdf",width = 10,height = 6)

prop_table_tot_test2 <- query_transfer@meta.data %>%
  # filter(robust_score_subclass != "uncertain") %>%
  # filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters_fix1,pathology_class,robust_score) %>%
  summarise(n = n()) %>%
  ungroup() %>% 
  group_by(seurat_clusters_fix1,pathology_class) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>% 
  ungroup()
# mutate(robust = case_when(robust_score!="uncertain"~"robust",
#                           T~robust_score))

prop_table_tot_test2 %>% 
  ggplot(aes(x=pathology_class,y=prop,col=seurat_clusters_fix1))+
  geom_point(position = position_jitter(width = 0.2))+
  facet_grid(~robust_score)+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45),strip.background = element_blank())
# ggsave("../../out/image/label_transfer_refWMImmune_queryBrainsphere_prop_table_tot_test2.pdf",width = 7,height = 3)

# heatmap for the assignament =============================================
prop_table_tot <- query_transfer@meta.data %>%
  # filter(robust_score_subclass != "uncertain") %>%
  # filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters_fix1,predicted.id) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters_fix1,predicted.id,prop) %>%
  pivot_wider(names_from = seurat_clusters_fix1,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.id") %>%
  as.matrix()

prop_table_robust <- query_transfer@meta.data %>%
  filter(robust_score != "uncertain") %>%
  group_by(seurat_clusters_fix1,predicted.id) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters_fix1,predicted.id,prop) %>%
  pivot_wider(names_from = seurat_clusters_fix1,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.id") %>%
  as.matrix()

# raw assignament
# pdf("../out/plot/129_heatmap_label_transfer_refWMImmune_queryIMMLYMsubset_raw_res0.9.pdf",height = 4,width = 5)
Heatmap(prop_table_tot,
        name = "prop", 
        column_title = "predicted WM id raw",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(11))
# dev.off()

# only on robust cell assignament
# pdf("../out/plot/129_heatmap_label_transfer_refWMImmune_queryIMMLYMsubset_robust_res0.9.pdf",height = 4,width = 5)
Heatmap(prop_table_robust,
        name = "prop", 
        column_title = "predicted WM id robust",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(11))
# dev.off()

# use the Jaccard score ==================================================
# similarly to the proportion of cells per annotation. I could also measure the cross cluster similarity per cell (how similar are the cluster from the query comapred to the annotation derived from the reference), using the Jaccard score.

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

# build the dataset for the correlatino plot cross the two annotation per cells
df_crossing <- crossing(id_ref = unique(query_transfer@meta.data$predicted.id),
                        id_query = unique(query_transfer@meta.data$seurat_clusters_fix1))

# df_crossing

# build the scatter plot
df_jaccard_score <- pmap(list(id_ref = df_crossing$id_ref,
                              id_query = df_crossing$id_query), function(id_ref,id_query){
                                
                                # calculate the jaccard score
                                a <- query_transfer@meta.data %>%
                                  rownames_to_column("barcode") %>%
                                  filter(predicted.id == id_ref) %>% pull(barcode)
                                
                                b <- query_transfer@meta.data %>%
                                  rownames_to_column("barcode") %>%
                                  filter(seurat_clusters_fix1 == id_query) %>% pull(barcode)
                                
                                jaccard_score <- jaccard(a,b)
                                
                                # build a data.frame
                                df <- data.frame("id_ref" = id_ref,
                                                 "id_query" = id_query,
                                                 "jaccard_score" = jaccard_score)
                                return(df)
                              }) %>%
  bind_rows()

# check the table
head(df_jaccard_score)

# shape it as a matrix
mat_jaccard_score <- df_jaccard_score %>%
  pivot_wider(names_from = id_ref,values_from = jaccard_score) %>%
  column_to_rownames("id_query")

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

# pdf("../out/plot/129_heatmap_jaccard_res0.9.pdf",height = 4,width = 5)
ht_02
# dev.off()

# Inverse label transfer: query over ref ----------------------------------

# Unimodal UMAP Projection
# In Seurat v4, we also enable projection of a query onto the reference UMAP structure. This can be achieved by computing the reference UMAP model and then calling MapQuery() instead of TransferData().
# this one run on the subset will overwrite the model

# ref_WM <- RunUMAP(ref_WM,
#                   dims = 1:30,
#                   reduction = "pca",
#                   return.model = TRUE)

query_transfer_inverse <- MapQuery(anchorset = test.anchors_0,
                                   reference = ref_WM,
                                   query = query_transfer,
                                   refdata = list(celltype = "cellID"),
                                   reference.reduction = "pca", reduction.model = "umap")

# wrangling ---------------------------------------------------------------
# add the coordinates to the metadata
query_transfer_inverse <- AddMetaData(query_transfer_inverse, metadata = query_transfer_inverse@reductions$ref.umap@cell.embeddings)

# save the object
saveRDS(query_transfer_inverse,"../out/object/129_query_transfer_QueryOverRef_res0.9.rds")


# EDA ---------------------------------------------------------------------

# plot the predicition score =============================================

# default plot 
# plot the raw prediction based on the reference dataset
DimPlot(query_transfer_inverse,label = T,group.by = "seurat_clusters_fix1",reduction = "ref.umap")

# plot the robust prediction
DimPlot(query_transfer_inverse,label = T,group.by = "robust_score",reduction = "ref.umap")

# costum plot
# average the position of the clusters
data2_transfer_avg_raw <- query_transfer_inverse@meta.data %>% group_by(seurat_clusters_fix1) %>% select(refUMAP_1, refUMAP_2) %>% summarize_all(mean)

# plot the raw prediction based on the reference dataset
# save the reference UMAPS
UMAP_ref_WM <- DimPlot(ref_WM,group.by = "seurat_clusters",label = T)
UMAP_ref_WM$data %>%
  ggplot(label= TRUE)+
  UMAP_ref_WM$layers[[1]]+
  UMAP_ref_WM$layers[[2]]+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5)))
# ggsave("../../out/image/129_UMAP_ref_WM.pdf",width = 5,height = 4)

# save the image with no label and no legend
UMAP_ref_WM$data %>%
  ggplot(label= F)+
  UMAP_ref_WM$layers[[1]]+
  # UMAP_ref_WM$layers[[2]]+
  theme_cowplot()+
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=5)))
# ggsave("../../out/image/129_UMAP_ref_WM_no_lab.pdf",width = 4,height = 4)

# plot all the cells in the query
ggplot(label= TRUE)+
  geom_point(data = ref_WM@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  geom_point(data = query_transfer_inverse@meta.data,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data2_transfer_avg_raw,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))
# ggsave("../out/plot/129_UMAP_label_transfer_refWMImmune_queryIMMLYMsubset_raw_inverse_res0.9.pdf",width = 6,height = 4)

# split by condition
ggplot(label= TRUE)+
  geom_point(data = ref_WM@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  geom_point(data = query_transfer_inverse@meta.data,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data2_transfer_avg_raw,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+facet_wrap(~pathology_class)+theme(strip.background = element_blank())
# ggsave("../out/plot/129_UMAP_label_transfer_refWMImmune_queryIMMLYMsubset_raw_inverse2_res0.9.pdf",width = 18,height = 12)

# save the image with no label and no legend
ggplot(label= F)+
  geom_point(data = ref_WM@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  geom_point(data = query_transfer_inverse@meta.data,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = data2_transfer_avg_raw,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+
  theme(legend.position = "none")
# ggsave("../out/plot/129_UMAP_label_transfer_refWMImmune_queryIMMLYMsubset_raw_inverse_no_lab_res0.9.pdf",width = 4,height = 4)

# explore the annotation
query_transfer_inverse@meta.data %>% 
  filter(robust_score != "uncertain") %>% 
  group_by(seurat_clusters_fix1,robust_score) %>% 
  summarise(n = n()) %>% 
  print(n=70)

# divide the dataset into uncertain and not based on the mapping score calculated bofore
data2_transfer_unc <- query_transfer_inverse@meta.data %>%
  filter(robust_score == "uncertain")
#
data2_transfer_defined <- query_transfer_inverse@meta.data %>%
  filter(robust_score != "uncertain")

# average the position of the clusters
data2_transfer_avg <- data2_transfer_defined %>% group_by(seurat_clusters_fix1) %>% select(refUMAP_1, refUMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE) +
  geom_point(data = ref_WM@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  # geom_point(data = data2_transfer_unc,aes(x = refUMAP_1,y = refUMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data2_transfer_defined,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data2_transfer_avg,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))
# ggsave("../../out/image/129_UMAP_label_transfer_refWMImmune_queryIMMLYMsubset_robust_inverse_res0.9.pdf",width = 8,height = 6)

# split by condition
ggplot(label= TRUE) +
  geom_point(data = ref_WM@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  # geom_point(data = data2_transfer_unc,aes(x = refUMAP_1,y = refUMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data2_transfer_defined,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data2_transfer_avg,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+facet_wrap(~pathology_class)+theme(strip.background = element_blank())
# ggsave("../out/plot/129_UMAP_label_transfer_refWMImmune_queryIMMLYMsubset_robust_inverse2_res0.9.pdf",width = 18,height = 12)

# print also without the labels
ggplot(label= F) +
  geom_point(data = ref_WM@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  # geom_point(data = data2_transfer_unc,aes(x = refUMAP_1,y = refUMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data2_transfer_defined,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = data2_transfer_avg,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+
  theme(legend.position = "none")
# ggsave("../out/plot/129_UMAP_label_transfer_refWMImmune_queryIMMLYMsubset_robust_inverse_no_lab_res0.9.pdf",width = 4,height = 4)

# df1 <- query_transfer_inverse@meta.data %>% 
#   group_by(pathology_class) %>% 
#   summarise(count_tot = n())
# 
# df2 <- query_transfer_inverse@meta.data %>% 
#   group_by(pathology_class) %>% 
#   summarise(count_robust = n())
# 
# left_join(df1,df2) %>% 
#   mutate(prop = count_robust/count_tot)

p1 <- ggplot(label= TRUE) +
  geom_point(data = ref_WM@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  # geom_point(data = data2_transfer_unc,aes(x = refUMAP_1,y = refUMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data2_transfer_defined,aes(x = refUMAP_1,y = refUMAP_2, col = predicted.celltype),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = data2_transfer_avg,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))

p2 <- UMAP_ref_WM$data %>%
  ggplot(label= F)+
  UMAP_ref_WM$layers[[1]]+
  # UMAP_ref_WM$layers[[2]]+
  theme_cowplot()+
  # theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=5)))

p1 + p2

# compare methods ---------------------------------------------------------
# check if the call from direct or inverse transfer produce the same annotation per cell
# notice that the inverse procedure is adding the embedding of the cells, and is also adding a predicted.celltype covariate.

# inverse transfer
query_transfer_inverse@meta.data$predicted.celltype

# direct transfer
query_transfer@meta.data$predicted.id

# comapre
table(query_transfer_inverse@meta.data$predicted.celltype,
      query_transfer@meta.data$predicted.id)

# the two covariates are the same.