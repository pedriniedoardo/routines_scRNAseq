# AIM ---------------------------------------------------------------------
# sample routine for the pseudobulk analysis.
# the following is a readaptation of different resources
# https://github.com/sib-swiss/single-cell-training/
# https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(GGally)
library(cowplot)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(DESeq2)
library(RNAseqQC)
library(limma)
library(ashr)
library(magick)
library(UpSetR)

# read in the data --------------------------------------------------------
# read in the sample dataset
data.combined <- readRDS("../out/object/100_ifnb_DonorStim.rds")

# add in one covariate the cell anntation and the stimulation status
data.combined$celltype.stim <- paste(data.combined$seurat_annotations, data.combined$stim, sep = "_")

# set the ident
Idents(data.combined) <- "celltype.stim"

# For this test I would focus on the subset only the presumed cells of interest cells for the test
scobj_subset <- subset(data.combined,subset = seurat_annotations %in% c("CD14 Mono"))
# DimPlot(scobj_subset_subset,label = T,raster = T,group.by = "stim")

# wrangling ---------------------------------------------------------------
# explore the full meta
scobj_subset@meta.data %>%
  group_by(donor_id_fix,stim) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = stim,values_from = n)

# aggregate the expression per sample per donor.id and stimulation
cts_sample_all <- AggregateExpression(object = scobj_subset,
                                      group.by = c("stim_donor"),
                                      assays = 'RNA',
                                      slot = "counts",
                                      return.seurat = FALSE)

# processing --------------------------------------------------------------
# extract
# 1. Get counts matrix
counts <- cts_sample_all$RNA %>%
  as.data.frame()

# 2. generate sample level metadata
LUT_sample <- scobj_subset@meta.data %>%
  group_by(stim_donor,donor_id,stim,seurat_annotations,celltype.stim) %>%
  summarise() %>%
  # harmonyze the naming of the donor_id
  mutate(stim_donor = str_replace(stim_donor,pattern = "_",replacement = "-"))

# match the LUT with the expression matrix
colData <- data.frame(sample.id = colnames(counts)) %>%
  # separate(samples,into = c(c("orig.ident","origin_01","pathology_class","origin_02")),remove = F,sep = "_") %>%
  left_join(LUT_sample,by = c("sample.id" = "stim_donor")) %>%
  # fix the donor_id
  mutate(donor_id = str_remove_all(donor_id,pattern = "-"))

# save matrix and metadata
saveRDS(counts,file = "../out/object/102_counts_mono_pBulk.rds")
saveRDS(colData,file = "../out/object/102_colData_mono_pBulk.rds")

# perform DESeq2 ----------------------------------------------------------

# build the model
treat <- colData$stim
block_donor <- colData$donor_id

design <- model.matrix(~ block_donor + treat)
# design <- model.matrix(~ treat)
colnames(design)[1] <- c("intercept")

# save the disign
saveRDS(design,"../out/object/102_design_mono_pBulk.rds")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = design)

# filter low aboundant features
feat_keep <- edgeR::filterByExpr(counts(dds), group = colData$stim)
dds_filter <- dds[feat_keep,]

# plot the raw number of reads per sample
colSums(counts(dds_filter)) %>%
  data.frame(tot_reads = .) %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x=sample,y=tot_reads)) + geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))

# scale the data
vds_filter <- vst(dds_filter, blind = F)

# clustering samples ------------------------------------------------------
# set seed to control random annotation colors
pdf("../out/plot/102_heatmap_SampleCluster_filterExp.pdf",width = 10,height = 6)
set.seed(1)
hm <- plot_sample_clustering(vds_filter,
                             anno_vars = c("sample.id","stim"),
                             distance = "euclidean")
draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

# PCA ---------------------------------------------------------------------
plot_vsd <- plotPCA(vds_filter,
                    intgroup = c("sample.id","stim")) +
  theme_bw()

plot_vsd$data %>%
  ggplot(aes(x=PC1,y=PC2,col=stim)) +
  geom_point(size =3,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd$labels[1]) + xlab(plot_vsd$labels[2])
ggsave("../out/plot/102_PCA_pseudobulk_filterExp.pdf",width = 6,height = 4)

# cell composition --------------------------------------------------------
# load the full dataset
scobj <- readRDS("../out/object/100_ifnb_DonorStim.rds")

# check the relative proportion of different cell types from different sample
sample_prop_wide <- scobj@meta.data %>%
  group_by(seurat_annotations,stim_donor,donor_id,stim) %>%
  summarise(n = n(),.groups = "drop") %>%
  # ensure all the combinations are presents
  complete(seurat_annotations,stim_donor, fill = list(n = 0)) %>%
  group_by(stim_donor) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  # make it wider
  dplyr::select(seurat_annotations,stim_donor,prop) %>%
  mutate(stim_donor = str_replace(stim_donor,pattern = "_",replacement = "-")) %>%
  pivot_wider(names_from = stim_donor,values_from = prop) %>%
  # pivot_wider(names_from = stim_donor,values_from = prop,values_fill = 0) %>%
  column_to_rownames("seurat_annotations") %>%
  # select only the samples in the dataset
  dplyr::select(colData$sample.id)

# confim the dimensions
colSums(sample_prop_wide)

# plot the data as heatmap
meta_sample_prop <- data.frame(colname = colnames(sample_prop_wide)) %>%
  left_join(colData,by=c("colname"="sample.id"))

column_meta_sample_prop <- HeatmapAnnotation(treat = meta_sample_prop$stim,
                                             col = list(treat = c("CTRL" = "blue",
                                                                  "STIM" = "yellow")))

ht2_shr <- Heatmap(sample_prop_wide, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                   name = "prop_cell_type",
                   column_title = "sample",
                   col = viridis::viridis(option = "turbo",n = 10),
                   
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_meta_sample_prop,show_row_names = T
                   # cluster_rows = F,
                   # right_annotation = row_ha,
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)
# pdf("../../out/image/heatmap_cellType_pseudobulk_ASTRO_ctrl_refCX.pdf",width = 7,height = 5)
draw(ht2_shr,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
# dev.off()

# scale the proportions by cell type
sample_prop_wide2 <- scobj@meta.data %>%
  group_by(seurat_annotations,stim_donor,donor_id,stim) %>%
  summarise(n = n(),.groups = "drop") %>%
  # ensure all the combinations are presents
  complete(seurat_annotations,stim_donor, fill = list(n = 0)) %>%
  group_by(stim_donor) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  ungroup() %>%
  group_by(seurat_annotations) %>%
  # mutate(test = (prop - mean(prop))/sd(prop)) %>%
  mutate(z_score = scale(prop)[,1]) %>%
  # make it wider
  dplyr::select(seurat_annotations,stim_donor,z_score) %>%
  mutate(stim_donor = str_replace(stim_donor,pattern = "_",replacement = "-")) %>%
  pivot_wider(names_from = stim_donor,values_from = z_score) %>%
  column_to_rownames("seurat_annotations") %>%
  # select only the samples in the dataset
  dplyr::select(colData$sample.id)

# confim the dimensions
rowSums(sample_prop_wide2)

# plot the data as heatmap
meta_sample_prop2 <- data.frame(colname = colnames(sample_prop_wide2)) %>%
  left_join(colData,by=c("colname"="sample.id"))

column_meta_sample_prop2 <- HeatmapAnnotation(treat = meta_sample_prop2$stim,
                                              col = list(treat = c("CTRL" = "blue",
                                                                   "STIM" = "yellow")))

ht2_shr2 <- Heatmap(sample_prop_wide2, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                    name = "zscore_cell_type",
                    column_title = "sample",
                    # col = viridis::viridis(option = "turbo",n = 10),
                    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                    
                    # row_names_gp = gpar(fontsize = 3),
                    top_annotation = column_meta_sample_prop2,show_row_names = T
                    # cluster_rows = F,
                    # right_annotation = row_ha,
                    # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)
# pdf("../../out/image/heatmap_cellType_pseudobulk_ASTRO_ctrl_refCX.pdf",width = 7,height = 5)
draw(ht2_shr2,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
# dev.off()

# run DE ------------------------------------------------------------------
# run DESeq2
ddsHTSeq_filter <- DESeq(dds_filter)

# Check the coefficients for the comparison
resultsNames(ddsHTSeq_filter)

# save the filtered object
saveRDS(ddsHTSeq_filter,"../out/object/102_ddsHTSeq_pseudobulk_filterExp.rds")

# print the contrast
resultsNames(ddsHTSeq_filter)

# save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast <- makeContrasts(CD14mono_STIM_vs_CTRL = treatSTIM,
                          levels = design)

# pull the results table
res_CD14mono_STIM <- results(ddsHTSeq_filter, contrast=contrast[,"CD14mono_STIM_vs_CTRL"],alpha = 0.05)
summary(res_CD14mono_STIM)

# add the gene symbols
df_res <-
  list(res_CD14mono_STIM = res_CD14mono_STIM) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsCTRL")
# save the table
df_res %>%
  write_tsv("../out/table/102_DE_pseudobulk_filterExp.tsv")

# shrink ------------------------------------------------------------------
res_CD14mono_STIM_shr <- lfcShrink(ddsHTSeq_filter, res = res_CD14mono_STIM, type = "ashr")

# save the table of DEGs
df_res_shr <-
  list(res_CD14mono_STIM_shr = res_CD14mono_STIM_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsCTRL")
# save the table
df_res_shr %>%
  write_tsv("../out/table/102_DE_pseudobulk_filterExp_shr.tsv")

# Another useful diagnostic plot is the histogram of the p values (figure below). This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
df_res %>%
  data.frame()%>%
  dplyr::filter(baseMean>1) %>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  facet_wrap(~conditionVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank())
ggsave("../out/plot/102_histogram_pvalue_pseudobulk_filterExp.pdf",width = 6,height = 4)

# volcano -----------------------------------------------------------------
# add the info of the genename
plot_volcano <- df_res %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
  #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
  #   size = 2,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  ggrepel::geom_text_repel(
    data = plot_volcano %>% 
      group_by(conditionVsCTRL) %>% 
      arrange(padj) %>% 
      dplyr::slice(1:30) %>% 
      dplyr::filter(abs(log2FoldChange)>1) %>% 
      dplyr::filter(padj<0.05),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../out/plot/102_vulcano_plot_pseudobulk_filterExp.pdf",width = 10,height = 10)

#
plot_volcano_shr <- df_res_shr %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano_shr %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = plot_volcano_shr %>% group_by(conditionVsCTRL) %>% arrange(padj) %>% dplyr::slice(1:10)%>% filter(abs(log2FoldChange)>1),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../out/plot/102_vulcano_plot_pseudobulk_filterExp_shr.pdf",width = 10,height = 10)

# MA plot -----------------------------------------------------------------
df_res %>%
  filter(!is.na(padj)) %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) +
  scale_x_log10() + scale_color_manual(values = c("gray","red")) +
  facet_wrap(~conditionVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")
ggsave("../out/plot/102_MA_plot_pseudobulk_filterExp.pdf",width = 10,height = 10)

# shr
df_res_shr %>%
  filter(!is.na(padj)) %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) +
  scale_x_log10() + scale_color_manual(values = c("gray","red")) +
  facet_wrap(~conditionVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")
ggsave("../out/plot/102_MA_plot_pseudobulk_filterExp_shr.pdf",width = 10,height = 10)

# heatmaps ----------------------------------------------------------------
# pull the scaled values
mat_filter <- assay(vds_filter) %>%
  as.data.frame() %>%
  # dplyr::select(contains(c("_0_","_6_"))) %>%
  as.matrix()

# the DEGs plot
DEG_2 <- df_res_shr %>%
  as.data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  dplyr::filter(col==1) %>%
  pull(symbol) %>%
  unique()

# mat_filter <- assay(vds_filter) %>%
#   data.frame() %>%
#   # dplyr::select(contains(c("_0_","_6_"))) %>%
#   as.matrix()

mat_shr <- mat_filter[rownames(vds_filter) %in% DEG_2, ]
mat2_shr <- (mat_shr - rowMeans(mat_shr))/rowSds(mat_shr,useNames = TRUE)
#
meta_sample <- data.frame(colname = colnames(mat2_shr)) %>%
  left_join(colData,by=c("colname"="sample.id"))

# make the column of the matrix more readable
colnames(mat2_shr) <- meta_sample$colname

column_ha_shr <- HeatmapAnnotation(treat = meta_sample$stim,
                                   # gender = meta_sample_ENDO$sex,
                                   # facility = meta_sample_ENDO$facility,
                                   col = list(treat = c("CTRL" = "blue",
                                                        "STIM" = "yellow")))

ht2_shr <- Heatmap(mat2_shr, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                   name = "exp",
                   column_title = "`test_shr",
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_ha_shr,show_row_names = F,
                   # cluster_rows = F,
                   # right_annotation = row_ha,
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)
pdf("../out/plot/102_heatmap_DEG_plot_pseudobulk_filterExp_shr.pdf",width = 7,height = 14)
draw(ht2_shr,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# # upset plot --------------------------------------------------------------
# # read in the table of DEGs
# df_res_shr <- read_tsv("../../out/table/DE_pseudobulk_ASTRO_ctrl_refCX_shr.tsv")
# 
# # build a list of common elements belonging to each set of fegs
# list_DE_up <- df_res_shr %>%
#   split(f = .$conditionVsCX) %>%
#   lapply(function(x){
#     x %>%
#       filter(padj < 0.05,log2FoldChange>1) %>%
#       pull(symbol) %>%
#       unique()
#   })
# 
# glimpse(list_DE_up)
# 
# list_DE_down <- df_res_shr %>%
#   split(f = .$conditionVsCX) %>%
#   lapply(function(x){
#     x %>%
#       filter(padj < 0.05,log2FoldChange<(-1)) %>%
#       pull(symbol) %>%
#       unique()
#   })
# glimpse(list_DE_down)
# 
# # try the upset plot version
# # library(UpSetR)
# pdf("../../out/image/upset_DEG_UP_plot_pseudobulk_ASTRO_refBASELINE_shr.pdf",width = 14,height = 7)
# upset(fromList(list_DE_up), order.by = "freq",nsets = 7)
# dev.off()
# 
# pdf("../../out/image/upset_DEG_DOWN_plot_pseudobulk_ASTRO_refBASELINE_shr.pdf",width = 14,height = 7)
# upset(fromList(list_DE_down), order.by = "freq",nsets = 7)
# dev.off()
# 
# # pull the intersections
# df1_UP <- lapply(list_DE_up,function(x){
#   data.frame(gene = x)
# }) %>%
#   bind_rows(.id = "path")
# 
# df1_DOWN <- lapply(list_DE_down,function(x){
#   data.frame(gene = x)
# }) %>%
#   bind_rows(.id = "path")
# 
# head(df1_UP)
# head(df1_DOWN)
# 
# df2_UP <- data.frame(gene=unique(unlist(list_DE_up)))
# df2_DOWN <- data.frame(gene=unique(unlist(list_DE_down)))
# 
# head(df2_UP)
# head(df2_DOWN)
# 
# df_int_UP <- lapply(df2_UP$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1_UP %>%
#     dplyr::filter(gene==x) %>%
#     arrange(path) %>%
#     pull("path") %>%
#     paste0(collapse = "|")
# 
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>%
#   bind_rows()
# 
# df_int_DOWN <- lapply(df2_DOWN$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1_DOWN %>%
#     dplyr::filter(gene==x) %>%
#     arrange(path) %>%
#     pull("path") %>%
#     paste0(collapse = "|")
# 
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>%
#   bind_rows()
# 
# df_int_UP %>%
#   write_tsv("../../out/table/upset_DEG_UP_plot_pseudobulk_ASTRO_refBASELINE_shr.tsv")
# 
# df_int_DOWN %>%
#   write_tsv("../../out/table/upset_DEG_DOWN_plot_pseudobulk_ASTRO_refBASELINE_shr.tsv")
# 
# head(df_int_UP,n=20)
# head(df_int_DOWN,n=20)
# 
# df_int_UP %>%
#   group_by(int) %>%
#   summarise(n=n()) %>%
#   arrange(desc(n))
# 
# df_int_DOWN %>%
#   group_by(int) %>%
#   summarise(n=n()) %>%
#   arrange(desc(n))

# PLOT DISPERSION ---------------------------------------------------------
pdf("../out/plot/102_dispersion_plot_pseudobulk_filterExp.pdf",width = 5,height = 5)
plotDispEsts(ddsHTSeq_filter)
dev.off()
