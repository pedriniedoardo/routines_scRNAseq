# AIM ---------------------------------------------------------------------
# This script is run to explore the dataset individually and define a consensus threshold for the mito and reads counts
# This is run on the post SoupX data correction.

# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# READ IN DATA ------------------------------------------------------------
folder <- "../../data/SopuX_out/"
id_sample <- dir(folder)

# data.frame(sample = id_sample) %>%
#   separate(col = "sample",into = c("treat","sample_id"),sep = "_",remove = F) %>%
#   write_csv("../data/LUT_samples.csv")

# load the LUT
LUT <- read_csv("../../data/LUT_sample.csv")

# do the preprocessing over all the dataset and save the objects
# x <- "test_neuron_auto"
list_datasc <- lapply(id_sample,function(x){
  # to track the processing of the progress of the lapply
  print(x)
  
  # read in the matrix
  data <- Read10X(data.dir = paste0(folder,x))
  
  # crete the object
  datasc <- CreateSeuratObject(counts = data, project = LUT %>%
                                 filter(sample_id == x) %>%
                                 pull(sample_name), min.cells = 20, min.features = 200)
  
  # add the metadata
  # datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
  # datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^mt-")
  datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")
  # add also the percentage of globin. in this dataset it is not meaningful as there is no blood
  # datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
  datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^Hb[^(p)]")
  
  # label the cells based on the mt reads content
  datasc$mt_bin <- datasc@meta.data %>%
    mutate(test = case_when(percent.mt < 1~"low",
                            percent.mt < 10~"mid",
                            T ~ "high")) %>%
    pull(test)
  
  # datasc$treat <- LUT %>%
  #   filter(sample_id == x) %>%
  #   pull(status)
  # 
  # datasc$ID <- LUT %>%
  #   filter(sample_id == x) %>%
  #   pull(sample_id)
  
  # add the filtering variable based on the fixed threshold. at this poit add some generic thresholds
  datasc$discard_threshold <- datasc@meta.data %>%
    # mutate(test = percent.mt > mito_thr | nFeature_RNA < featureLow_thr | nFeature_RNA > featureHigh_thr) %>% 
    mutate(test = percent.mt > 15 | nFeature_RNA < 1000 | nFeature_RNA > 6000) %>% 
    pull(test)
  
  # add the filtering variable based on the adaptive threshold, multivariable
  # library(robustbase)
  # library(scater)
  stats <- cbind(log10(datasc@meta.data$nCount_RNA),
                 log10(datasc@meta.data$nFeature_RNA),
                 datasc@meta.data$percent.mt)
  
  outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
  multi.outlier <- isOutlier(outlying, type = "higher")
  datasc$discard_multi <- as.vector(multi.outlier)
  
  # add the filtering variable based on the adaptive threshold single variable
  high_QC_mito <- isOutlier(datasc@meta.data$percent.mt, type="high", log=TRUE)
  QC_features <- isOutlier(datasc@meta.data$nFeature_RNA, type="both", log=TRUE)
  
  datasc$discard_single <- high_QC_mito | QC_features
  
  return(datasc)
}) %>%
  setNames(id_sample)

# confirm the class of the objects ----------------------------------------
lapply(list_datasc, function(x){
  class(x@assays$RNA)
})

lapply(list_datasc, function(x){
  dim(x)
})

# plot QC -----------------------------------------------------------------
# extract the metadata from each dataset
meta_total <- lapply(list_datasc, function(x){
  x@meta.data %>%
    rownames_to_column("barcode") %>%
    mutate(barcode = paste0(barcode,"|",orig.ident)) %>%
    mutate(tot_nFeature = rownames(x) %>% length())
}) %>%
  bind_rows(.id = "dataset")

meta_total %>%
  write_tsv("../../out/R_analysis/table/01_meta_datasc_beforeQC_V5.tsv")

# meta_total <-read_tsv("../../out/table/meta_datasc_fix_filter_norm_total.tsv")

# how many cells are considered outliers
meta_total %>%
  dplyr::count(dataset,discard_threshold)

meta_total %>%
  dplyr::count(dataset,discard_multi)

meta_total %>%
  dplyr::count(dataset,discard_single)

# fixed threshold scatter nFeature vs percent.mt
meta_total %>%
  ggplot(aes(y = percent.mt,x = nFeature_RNA,col=discard_threshold)) + geom_point(alpha=0.3) +
  facet_wrap(~dataset) +
  theme_bw() +
  theme(strip.background = element_blank())
# save the plot
ggsave("../../out/R_analysis/plot/01_fixed_scatter_feature_mito_V5.pdf",width = 8,height = 4)

meta_total %>%
  ggplot(aes(y = percent.mt,x = nFeature_RNA,col=discard_multi)) +
  geom_point(alpha=0.3) +
  facet_wrap(~dataset) +
  theme_bw() +
  theme(strip.background = element_blank())
# save the plot
ggsave("../../out/R_analysis/plot/01_adaptiveMulti_scatter_feature_mito_V5.pdf",width = 8,height = 4)

meta_total %>%
  ggplot(aes(y = percent.mt,x = nFeature_RNA,col=discard_single)) +
  geom_point(alpha=0.3) +
  facet_wrap(~dataset) +
  theme_bw() +
  theme(strip.background = element_blank())
# save the plot
ggsave("../../out/R_analysis/plot/01_adaptiveSingle_scatter_feature_mito_V5.pdf",width = 8,height = 4)

#
meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  ggplot(aes(x=factor(dataset),y=value)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.01) +
  facet_wrap(~var,scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  theme(strip.background = element_blank()) +
  scale_y_continuous(trans = "log1p")
# save the plot
ggsave("../../out/R_analysis/plot/01_fixed_boxplot_reads_V5.pdf",width = 6,height = 6)

#
meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "percent.mt") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
  facet_wrap(dataset~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(0,5,10,20,40,100)) +
  geom_vline(xintercept = c(10),col="red",linetype="dashed") +
  annotate("rect", xmin=0, xmax=10, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())
# save the plot
ggsave("../../out/R_analysis/plot/01_fixed_histo_mito_V5.pdf",width = 8,height = 4)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nFeature_RNA") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(dataset~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(500,1000,2000,4000,8000,16000)) +
  geom_vline(xintercept = c(1000,7000),col="red",linetype="dashed") +
  annotate("rect", xmin=1000, xmax=7000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())
# save the plot
ggsave("../../out/R_analysis/plot/01_fixed_histo_features_V5.pdf",width = 8,height = 4)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nFeature_RNA") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(dataset~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p") +
  geom_vline(xintercept = c(1000,5000),col="red",linetype="dashed") +
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nCount_RNA") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(dataset~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p") +
  # geom_vline(xintercept = c(500,5000),col="red",linetype="dashed") +
  # annotate("rect", xmin=500, xmax=5000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())
# save the plot
ggsave("../../out/R_analysis/plot/01_fixed_histo_counts_V5.pdf",width = 8,height = 4)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.ribo,percent.globin)) %>%
  filter(var == "percent.ribo") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(dataset~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p") +
  # geom_vline(xintercept = c(500,5000),col="red",linetype="dashed") +
  # annotate("rect", xmin=500, xmax=5000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())
# save the plot
ggsave("../../out/R_analysis/plot/01_fixed_histo_ribo_V5.pdf",width = 8,height = 4)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.ribo,percent.globin)) %>%
  filter(var == "percent.globin") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(dataset~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p") +
  # geom_vline(xintercept = c(500,5000),col="red",linetype="dashed") +
  # annotate("rect", xmin=500, xmax=5000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())
# save the plot
ggsave("../../out/R_analysis/plot/01_fixed_histo_globin_V5.pdf",width = 8,height = 4)

# color the bins for the amount of reads
meta_total %>%
  ggplot(aes(x = nCount_RNA,y = nFeature_RNA,col=mt_bin)) + geom_point(alpha=0.3) + facet_grid(orig.ident~mt_bin,scales = "free_y")+theme_bw() +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(strip.background = element_blank())
# save the plot
# ggsave("out/image/fixed_threshold/fixed_scatter_mito.pdf",width = 10,height = 6)