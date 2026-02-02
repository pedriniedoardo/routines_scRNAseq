# AIM ---------------------------------------------------------------------
# run a statistical test on the proportion differences across clusters, between conditions

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(cacoa)
library(conos)
library(sccore)
library(coda.base)
library(psych)

# build the cacoa object from seurat --------------------------------------
# load the seurat object
so_ref <- readRDS("../../out/object/sobj_total_full_h.rds")
so_test <- readRDS("../../out/object/sobj_total_remove_h.rds")

DimPlot(so_ref,label = T,group.by = "test_cellid",raster=T,split.by = "test_disease")+plot_annotation("so_ref")
DimPlot(so_test,label = T,group.by = "test_cellid",raster=T,split.by = "test_disease")+plot_annotation("so_test")

# wrangling ---------------------------------------------------------------
# extract the meta
meta_ref <- so_ref@meta.data %>%
  rownames_to_column()

meta_test <- so_test@meta.data %>%
  rownames_to_column()

# check dimensions --------------------------------------------------------
# confirm the numbers from the reference dataset
head(meta_ref)
dim(meta_ref)

head(meta_test)
dim(meta_test)

# run the proportion test diagnosis ---------------------------------------
# Run propeller testing for cell type proportion differences between the groups.
# cluster is the cluster/celltype id
# sample is the reference id of the biological replicates
# group is the grouping id

# ref dataset
properller_out_ref <- propeller(clusters = meta_ref$test_cellid,
                                sample = meta_ref$test_donor,
                                group = meta_ref$test_disease)

properller_out_ref %>%
  rownames_to_column("cell_id") %>%
  write_tsv("../../out/table/02_propeller_out_ref_seuratObject_simulatedData.tsv")

# test dataset
properller_out_test <- propeller(clusters = meta_test$test_cellid,
                                 sample = meta_test$test_donor,
                                 group = meta_test$test_disease)

properller_out_test %>%
  rownames_to_column("cell_id") %>%
  write_tsv("../../out/table/02_propeller_out_test_seuratObject_simulatedData.tsv")

# plotting diagnosis ------------------------------------------------------
# reference
df_summary_ref <- meta_ref %>% 
  group_by(test_cellid,
           test_donor,
           test_disease) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(test_donor) %>% 
  mutate(tot = sum(n),
         prop = n/tot) %>%
  ungroup()

# test
df_summary_test <- meta_test %>% 
  group_by(test_cellid,
           test_donor,
           test_disease) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(test_donor) %>% 
  mutate(tot = sum(n),
         prop = n/tot) %>%
  ungroup()

# plot 01
df_summary_ref %>%
  ggplot(aes(x=test_disease,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~test_cellid,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))+
  ggtitle("summary so_ref")
# ggsave("../../out/plot/manualClean/propeller_plot01_ref.pdf",width = 10,height = 10)

df_summary_test %>%
  ggplot(aes(x=test_disease,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~test_cellid,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))+
  ggtitle("summary so_test")
# ggsave("../../out/plot/manualClean/propeller_plot01_test.pdf",width = 10,height = 10)

# plot 02
df_summary_ref %>%
  ggplot() +
  geom_boxplot(aes(x=test_cellid,y=prop,color=test_disease),outlier.shape = NA) +
  geom_point(aes(x=test_cellid,y=prop,color=test_disease),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
# ggsave("../../out/plot/manualClean/propeller_plot02_diagnosis_cellid.pdf",width = 8,height = 5)

df_summary_test %>%
  ggplot() +
  geom_boxplot(aes(x=test_cellid,y=prop,color=test_disease),outlier.shape = NA) +
  geom_point(aes(x=test_cellid,y=prop,color=test_disease),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
# ggsave("../../out/plot/manualClean/propeller_plot02_test.pdf",width = 8,height = 5)

# plot 03
df_summary_ref %>%
  group_by(test_disease,test_cellid) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(test_disease) %>%
  mutate(tot = sum(n),
         prop = n/tot) %>%
  ggplot() +
  geom_col(aes(x=test_disease,y=prop,fill=test_cellid))+
  theme_cowplot()

df_summary_test %>%
  group_by(test_disease,test_cellid) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(test_disease) %>%
  mutate(tot = sum(n),
         prop = n/tot) %>%
  ggplot() +
  geom_col(aes(x=test_disease,y=prop,fill=test_cellid))+
  theme_cowplot()
