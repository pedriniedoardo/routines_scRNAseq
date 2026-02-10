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
options(Seurat.object.assay.version = "v3")

# READ IN DATA ------------------------------------------------------------
id_sample <- dir("../data/test_introns/SoupX_default_cellranger7/")

# load the LUT
LUT <- read_csv("../data/test_introns/LUT_samples.csv")

# do the preprocessing over all the dataset and save the objects
# x <- "sample_untreated_Wintron"
list_datasc <- lapply(id_sample,function(x){
  # to track the processing of the progress of the lapply
  print(x)
  
  # read in the matrix
  data <- Read10X(data.dir = paste0("../data/test_introns/SoupX_default_cellranger7/",x))
  
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
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
  datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  # add also the percentage of globin. in this dataset it is not meaningful as there is no blood
  datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
  
  # label the cells based on the mt reads content
  datasc$mt_bin <- datasc@meta.data %>%
    mutate(test = case_when(percent.mt < 1~"low",
                            percent.mt < 10~"mid",
                            T ~ "high")) %>%
    pull(test)
  
  datasc$treat <- LUT %>%
    filter(sample == x) %>%
    pull(treat)
  
  datasc$test_intron <- LUT %>%
    filter(sample == x) %>%
    pull(test_intron)
  
  datasc$doxy <- LUT %>%
    filter(sample == x) %>%
    pull(doxy)
  
  datasc$exposure <- LUT %>%
    filter(sample == x) %>%
    pull(exposure)
  
  datasc$ID <- LUT %>%
    filter(sample == x) %>%
    pull(ID)
  
  # add the filtering variable based on the fixed threshold
  datasc$test <- datasc@meta.data %>%
    mutate(test = percent.mt < 15 & nFeature_RNA > 1000 & nFeature_RNA < 6000) %>% 
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

# plot QC -----------------------------------------------------------------
# extract the metadata from each dataset
meta_total <- lapply(list_datasc, function(x){
  x@meta.data %>%
    rownames_to_column("barcode") %>%
    mutate(barcode = paste0(barcode,"|",orig.ident))
}) %>%
  bind_rows(.id = "dataset")

meta_total %>%
  write_tsv("../out/test_introns/table/meta_datasc_beforeQC_V3.tsv")

# meta_total <-read_tsv("../../out/table/meta_datasc_fix_filter_norm_total.tsv")

# how many cells are considered outliers
meta_total %>%
  dplyr::count(ID,test)

meta_total %>%
  dplyr::count(ID,not_outlier)

# fixed threshold scatter nFeature vs percent.mt
meta_total %>%
  ggplot(aes(y = percent.mt,x = nFeature_RNA,col=test)) + geom_point(alpha=0.3) +
  facet_wrap(~ID) +
  theme_bw() +
  theme(strip.background = element_blank())
# save the plot
ggsave("../out/test_introns/plot/fixed_scatter_feature_mito_V3.pdf",width = 8,height = 4)

meta_total %>%
  ggplot(aes(y = percent.mt,x = nFeature_RNA,col=not_outlier)) +
  geom_point(alpha=0.3) +
  facet_wrap(~ID) +
  theme_bw() +
  theme(strip.background = element_blank())
# save the plot
ggsave("../out/test_introns/plot/adaptive_scatter_feature_mito_V3.pdf",width = 8,height = 4)

#
meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  ggplot(aes(x=ID,y=value)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.01) +
  facet_wrap(~var,scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  theme(strip.background = element_blank())
# save the plot
ggsave("../out/test_introns/plot/fixed_boxplot_reads_V3.pdf",width = 6,height = 6)

#
meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "percent.mt") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
  facet_wrap(ID~var,scales = "free") +
  theme_bw() +
  scale_x_log10() +
  geom_vline(xintercept = c(15),col="red",linetype="dashed") +
  annotate("rect", xmin=0, xmax=15, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())
# save the plot
ggsave("../out/test_introns/plot/fixed_histo_mito_V3.pdf",width = 8,height = 4)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nFeature_RNA") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(ID~var,scales = "free") +
  theme_bw() +
  scale_x_log10() +
  geom_vline(xintercept = c(1000,6000),col="red",linetype="dashed") +
  annotate("rect", xmin=1000, xmax=6000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())
# save the plot
ggsave("../out/test_introns/plot/fixed_histo_features_V3.pdf",width = 8,height = 4)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nCount_RNA") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(ID~var,scales = "free") +
  theme_bw() +
  scale_x_log10() +
  # geom_vline(xintercept = c(500,5000),col="red",linetype="dashed") +
  # annotate("rect", xmin=500, xmax=5000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())
# save the plot
ggsave("../out/test_introns/plot/fixed_histo_counts_V3.pdf",width = 8,height = 4)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.ribo,percent.globin)) %>%
  filter(var == "percent.ribo") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(ID~var,scales = "free") +
  theme_bw() +
  scale_x_log10() +
  # geom_vline(xintercept = c(500,5000),col="red",linetype="dashed") +
  # annotate("rect", xmin=500, xmax=5000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())
# save the plot
ggsave("../out/test_introns/plot/fixed_histo_ribo_V3.pdf",width = 8,height = 4)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.ribo,percent.globin)) %>%
  filter(var == "percent.globin") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(ID~var,scales = "free") +
  theme_bw() +
  scale_x_log10() +
  # geom_vline(xintercept = c(500,5000),col="red",linetype="dashed") +
  # annotate("rect", xmin=500, xmax=5000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())
# save the plot
ggsave("../out/test_introns/plot/fixed_histo_globin_V3.pdf",width = 8,height = 4)

# color the bins for the amount of reads
meta_total %>%
  ggplot(aes(x = nCount_RNA,y = nFeature_RNA,col=mt_bin)) + geom_point(alpha=0.3) + facet_grid(orig.ident~mt_bin,scales = "free_y")+theme_bw() +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(strip.background = element_blank())
# save the plot
# ggsave("out/image/fixed_threshold/fixed_scatter_mito.pdf",width = 10,height = 6)