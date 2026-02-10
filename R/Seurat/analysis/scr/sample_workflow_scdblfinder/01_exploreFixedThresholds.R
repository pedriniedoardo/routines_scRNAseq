# AIM ---------------------------------------------------------------------
# This script is run to explore the dataset individually and define a consensus threshold for the mito and reads counts
# This is run on the post SoupX data correction.

# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(ggrastr)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# READ IN DATA ------------------------------------------------------------
folder <- "data/raw/cellranger8/SoupX/"
id_sample <- dir(folder)

# data.frame(sample = id_sample) %>%
#   separate(col = "sample",into = c("treat","sample_id"),sep = "_",remove = F) %>%
#   write_csv("../data/LUT_samples.csv")

# load the LUT
LUT <- read_csv("data/LUT_samples.csv")

# create a proposed input for the QC filters
LUT_QC_before <- data.frame(sample_name = id_sample,
           featureLow_thr = 500,
           featureHigh_thr = 6000,
           mito_thr = 10)

# save the table as reference
LUT_QC_before %>%
  write_csv("data/LUT_QC_before.csv")

# save the table to be update afterwards
# LUT_QC_before %>%
#   write_csv("data/LUT_QC_after.csv")

# do the preprocessing over all the dataset and save the objects
# x <- "s1"
list_datasc <- lapply(id_sample,function(x){
  # to track the processing of the progress of the lapply
  print(x)
  
  # read in the matrix
  data <- Read10X(data.dir = paste0(folder,x))
  
  # crete the object
  datasc <- CreateSeuratObject(counts = data, project = LUT %>%
                                 filter(sample == x) %>%
                                 pull(sample), 
                               # remove low expressed features
                               # in this case keep all the features before the integration
                               min.cells = 0,
                               # remove low content barcodes
                               min.features = 200)
  
  # add the metadata
  # datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
  # datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
  datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  # add also the percentage of globin. in this dataset it is not meaningful as there is no blood
  # datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
  datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
  
  # label the cells based on the mt reads content
  datasc$mt_bin <- datasc@meta.data %>%
    mutate(test = case_when(percent.mt < 1~"low",
                            percent.mt < 10~"mid",
                            T ~ "high")) %>%
    pull(test)
  
  datasc$pathology <- LUT %>%
    filter(sample == x) %>%
    pull(pathology)
  
  datasc$location <- LUT %>%
    filter(sample == x) %>%
    pull(location)
  
  datasc$tissue <- LUT %>%
    filter(sample == x) %>%
    pull(tissue)
  
  # add the filtering variable based on the fixed threshold
  # use the proposed filters from the LUT_QC_before
  df_QC <- LUT_QC_before %>%
    filter(sample_name == x)
  
  # add the filters
  datasc$test <- datasc@meta.data %>%
    mutate(test = percent.mt < df_QC$mito_thr & nFeature_RNA > df_QC$featureLow_thr & nFeature_RNA < df_QC$featureHigh_thr) %>% 
    pull(test)
  
  # add the filtering variable based on the adaptive threshold
  # library(robustbase)
  # library(scater)
  stats <- cbind(log10(datasc@meta.data$nCount_RNA),
                 log10(datasc@meta.data$nFeature_RNA),
                 datasc@meta.data$percent.mt)
  
  outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
  multi.outlier <- isOutlier(outlying, type = "higher")
  
  datasc$not_outlier <- !as.vector(multi.outlier)
  
  return(datasc)
}) %>%
  setNames(id_sample)

# confirm the class of the objects ----------------------------------------
lapply(list_datasc, function(x){
  class(x@assays$RNA)
})

# plot QC before filtering ------------------------------------------------
# extract the metadata from each dataset
meta_total <- lapply(list_datasc, function(x){
  x@meta.data %>%
    rownames_to_column("barcode") %>%
    mutate(barcode = paste0(barcode,"|",orig.ident))
}) %>%
  bind_rows(.id = "dataset") %>%
  left_join(y = LUT_QC_before,by = c("dataset" = "sample_name"))

meta_total %>%
  write_tsv("out/table/analysis_R45/meta_datasc_beforeQC_V5.tsv")

# confirm the addition of the filters
meta_total %>%
  group_by(dataset,featureLow_thr,featureHigh_thr,mito_thr) %>%
  summarise(n = n())

# how many cells are considered outliers
meta_total %>%
  dplyr::count(orig.ident,test)

meta_total %>%
  dplyr::count(orig.ident,not_outlier)

# fixed threshold scatter nFeature vs percent.mt
p01 <- meta_total %>%
  ggplot(aes(y = percent.mt,x = nFeature_RNA,col=test)) + geom_point(alpha=0.3) +
  facet_wrap(~orig.ident) +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# save the plot
ggsave(plot = p01,
       filename = "out/plot/analysis_R45/fixed_scatter_feature_mito_beforeQC_V5.pdf",width = 8,height = 4)

# adaptive threshold scatter nFeature vs percent.mt
p02 <- meta_total %>%
  ggplot(aes(y = percent.mt,x = nFeature_RNA,col=not_outlier)) +
  geom_point(alpha=0.3) +
  facet_wrap(~orig.ident) +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# save the plot
ggsave(plot = p02,
       filename = "out/plot/analysis_R45/adaptive_scatter_feature_mito_beforeQC_V5.pdf",width = 8,height = 4)

# plot the QC parameters as boxplots
p03 <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  ggplot(aes(x=factor(orig.ident),y=value)) +
  geom_violin() +
  # vectorial solution might be too high
  # geom_jitter(width = 0.2, alpha = 0.01) +
  # use raster solution to reduce the size
  geom_jitter_rast(width = 0.2, alpha = 0.01, raster.dpi = 300) +
  facet_wrap(~var,scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  theme(strip.background = element_blank()) +
  scale_y_continuous(trans = "log1p")
# save the plot
ggsave(plot = p03,
       filename = "out/plot/analysis_R45/fixed_boxplot_reads_beforeQC_V5.pdf",width = 6,height = 6)

# plot the distribution of the mito porp reads with cut-off
p04 <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "percent.mt") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
  facet_wrap(orig.ident~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p") +
  # 1. Use geom_rect instead of annotate
  geom_rect(
    data = LUT_QC_before %>% dplyr::rename(orig.ident = sample_name),
    aes(xmin = 0, xmax = mito_thr, ymin = 0, ymax = Inf, x = NULL),
    alpha = 0.1, 
    fill = "red") +
  # 2. Use the same data for the dashed lines
  geom_vline(
    data = LUT_QC_before %>% dplyr::rename(orig.ident = sample_name),
    aes(xintercept = mito_thr), 
    col = "red", linetype = "dashed"
  ) +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# save the plot
ggsave(plot = p04,
       filename = "out/plot/analysis_R45/fixed_histo_mito_beforeQC_V5.pdf",width = 8,height = 4)

# plot the distribution of the nfeature with cut-off
p05 <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nFeature_RNA") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(ID~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p") +
  # 1. Use geom_rect instead of annotate
  geom_rect(
    data = LUT_QC_before %>% dplyr::rename(orig.ident = sample_name),
    aes(xmin = featureLow_thr,xmax = featureHigh_thr,ymin = 0, ymax = Inf,x = NULL),
    alpha = 0.1,
    fill = "red") +
  # 2. Use the same data for the dashed lines
  geom_vline(
    data = LUT_QC_before %>% dplyr::rename(orig.ident = sample_name),
    aes(xintercept = featureLow_thr), 
    col = "red", linetype = "dashed"
  ) +
  geom_vline(
    data = LUT_QC_before %>% dplyr::rename(orig.ident = sample_name),
    aes(xintercept = featureHigh_thr), 
    col = "red", linetype = "dashed"
  ) +
  facet_wrap(orig.ident ~ var, scales = "free") +
  scale_x_continuous(trans = "log1p") +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# save the plot
ggsave(plot = p05,
       filename = "out/plot/analysis_R45/fixed_histo_features_beforeQC_V5.pdf",width = 8,height = 4)

# plot the distribution of the counts
p06 <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nCount_RNA") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(orig.ident~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p") +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# save the plot
ggsave(plot = p06,
       filename = "out/plot/analysis_R45/fixed_histo_counts_beforeQC_V5.pdf",width = 8,height = 4)

# plot the distribution of the ribo prop reads
p07 <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.ribo,percent.globin)) %>%
  filter(var == "percent.ribo") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(orig.ident~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p") +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# save the plot
ggsave(plot = p07,
       filename = "out/plot/analysis_R45/fixed_histo_ribo_beforeQC_V5.pdf",width = 8,height = 4)

# plot the distribution of the globin prop reads
p08 <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.ribo,percent.globin)) %>%
  filter(var == "percent.globin") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(orig.ident~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p") +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# save the plot
ggsave(plot = p08,
       filename = "out/plot/analysis_R45/fixed_histo_globin_beforeQC_V5.pdf",width = 8,height = 4)

# color the bins for the amount of reads
meta_total %>%
  ggplot(aes(x = nCount_RNA,y = nFeature_RNA,col=mt_bin)) + geom_point(alpha=0.3) + facet_grid(orig.ident~mt_bin,scales = "free_y")+theme_bw() +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(strip.background = element_blank())
# save the plot
# ggsave("out/image/fixed_threshold/fixed_scatter_mito.pdf",width = 10,height = 6)

# plot QC after filtering -------------------------------------------------
# produce the plots with filters after the QC evaluation
# evantually modify the values in the QC
LUT_QC_after <- read_csv("data/LUT_QC_after.csv")

# plot the distribution of the mito porp reads with cut-off
p09 <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "percent.mt") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
  facet_wrap(orig.ident~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p") +
  # 1. Use geom_rect instead of annotate
  geom_rect(
    data = LUT_QC_before %>% dplyr::rename(orig.ident = sample_name),
    aes(xmin = 0, xmax = mito_thr, ymin = 0, ymax = Inf, x = NULL),
    alpha = 0.1, 
    fill = "red") +
  # 2. Use the same data for the dashed lines
  geom_vline(
    data = LUT_QC_before %>% dplyr::rename(orig.ident = sample_name),
    aes(xintercept = mito_thr), 
    col = "red", linetype = "dashed"
  ) +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# save the plot
ggsave(plot = p09,
       filename = "out/plot/analysis_R45/fixed_histo_mito_afterQC_V5.pdf",width = 8,height = 4)

# plot the distribution of the nfeature with cut-off
p10 <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nFeature_RNA") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(ID~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p") +
  # 1. Use geom_rect instead of annotate
  geom_rect(
    data = LUT_QC_after %>% dplyr::rename(orig.ident = sample_name),
    aes(xmin = featureLow_thr,xmax = featureHigh_thr,ymin = 0, ymax = Inf,x = NULL),
    alpha = 0.1,
    fill = "red") +
  # 2. Use the same data for the dashed lines
  geom_vline(
    data = LUT_QC_after %>% dplyr::rename(orig.ident = sample_name),
    aes(xintercept = featureLow_thr), 
    col = "red", linetype = "dashed"
  ) +
  geom_vline(
    data = LUT_QC_after %>% dplyr::rename(orig.ident = sample_name),
    aes(xintercept = featureHigh_thr), 
    col = "red", linetype = "dashed"
  ) +
  facet_wrap(orig.ident ~ var, scales = "free") +
  scale_x_continuous(trans = "log1p") +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# save the plot
ggsave(plot = p10,
       filename = "out/plot/analysis_R45/fixed_histo_features_afterQC_V5.pdf",width = 8,height = 4)
