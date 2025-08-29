# AIM ---------------------------------------------------------------------
# sample plotting of a set of specific genes

# libraries ---------------------------------------------------------------
library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)
library(tidyverse)
library(ggrepel)
library(scales)
library(RColorBrewer)
library(SeuratWrappers)
library(dittoSeq)
library(clustree)
library(pals)
library(patchwork)
library(magick)
library(homologene)
library(SeuratWrappers)
library(presto)
library(GGally)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
# options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
data.combined <- readRDS("../out/object/100_ifnb_DonorStim.rds")
Idents(data.combined) <- "seurat_annotations"

GOI <- c("CCL3","CD8A")

# wrangling ---------------------------------------------------------------
# tidy up the metadata
# data.combined$donor_id_fix <- str_replace_all(data.combined$donor_id_fix,pattern = "-","_")
data.combined$seurat_annotations <- str_replace_all(data.combined$seurat_annotations,pattern = "\\s","-")

# EDA ---------------------------------------------------------------------
# check the global expression
VlnPlot(object = data.combined,features = GOI,group.by = "seurat_annotations")
# ggsave("../../out/plot/VlnProcr_cluster.pdf",width = 6,height = 4)

# manual plotting ---------------------------------------------------------
# extract the expression data
df_exp <- FetchData(data.combined, vars = GOI,layer = "data") |> 
  rownames_to_column("barcodes") |> 
  pivot_longer(names_to = "gene",values_to = "count",-barcodes) |> 
  # try to min/max normalize the count varaible per gene in order to rescale the difference in terms of expression
  group_by(gene) %>%
  # threshold of positiveness is based on the distriubtion of the expression of the signal in tihs case
  mutate(norm_min_max = ((count - min(count))/(max(count)-min(count))),
         exp_cat = case_when(count > 0~"pos",
                             T~"neg")) %>%
  ungroup() %>%
  mutate(count_fix = count + rnorm(nrow(.))/100000)
  # separate(barcodes,into = c("barcode","barcode_id"),sep = "-",remove = F)

head(df_exp)

exp_sc_wide <- df_exp %>%
  select(barcodes,gene,count_fix) %>%
  pivot_wider(names_from = gene,values_from = count_fix)

ggpairs(exp_sc_wide[,-1]) + theme_bw()

# average expression
# build the grouping variable
group_id_treat <- paste(data.combined@meta.data$donor_id_fix,
                        data.combined@meta.data$seurat_annotations,
                        data.combined@meta.data$stim,
                        sep = "|")
# add it to the meta
data.combined$group_id_treat <- group_id_treat

# set the idents
Idents(data.combined) <- "group_id_treat"
# DefaultAssay(data.combined) <- "RNA"

df_average_treat <- AverageExpression(data.combined,features = GOI)$RNA %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group_id",values_to = "avg_exp",-gene)

# build the pattern to extract the metadata
pattern_cellid <- paste0(unique(data.combined@meta.data$seurat_annotations),collapse = "|")
pattern_sample <- paste0(unique(data.combined@meta.data$donor_id_fix),collapse = "|")
pattern_treat <- paste0(unique(data.combined@meta.data$stim),collapse = "|")

df_average_treat %>%
  mutate(cellid = str_extract_all(group_id,pattern = pattern_cellid) %>% unlist()) %>%
  mutate(sample = str_extract_all(group_id,pattern = pattern_sample) %>% unlist()) %>%
  mutate(treat = str_extract_all(group_id,pattern = pattern_treat) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = cellid,y=avg_exp)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),shape = 1)+
  theme_bw() +
  facet_wrap(treat~gene)+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))
ggsave("../out/plot/104_AvgExpGOI_treat.pdf",width = 8,height = 8)

# try to color by treat rather than splitting
df_average_treat %>%
  mutate(cellid = str_extract_all(group_id,pattern = pattern_cellid) %>% unlist()) %>%
  mutate(sample = str_extract_all(group_id,pattern = pattern_sample) %>% unlist()) %>%
  mutate(treat = str_extract_all(group_id,pattern = pattern_treat) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = cellid,y=avg_exp,col=treat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),shape = 1)+
  theme_bw() +
  facet_wrap(~gene)+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))
ggsave("../out/plot/104_AvgExpGOI_treat_02.pdf",width = 8,height = 4)

# test if it is significant
df_average_treat %>%
  mutate(cellid = str_extract_all(group_id,pattern = pattern_cellid) %>% unlist()) %>%
  mutate(sample = str_extract_all(group_id,pattern = pattern_sample) %>% unlist()) %>%
  mutate(treat = str_extract_all(group_id,pattern = pattern_treat) %>% unlist()) %>%
  lm(data = .,avg_exp~cellid+treat) %>%
  summary()
