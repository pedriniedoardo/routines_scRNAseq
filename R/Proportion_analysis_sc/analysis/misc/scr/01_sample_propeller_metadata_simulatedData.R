# AIM ---------------------------------------------------------------------
# run a statistical test on the proportion differences across clusters, between conditions

# libraries ---------------------------------------------------------------
library(tidyverse)
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)
library(finalfit)

# read in the data --------------------------------------------------------
# read in the summary simualted metadata
df <- read_tsv("../../out/table/sample_metadata_simulatedData.tsv")

# generate a fake metadata based on the summary count
# ref data
meta_ref <- df %>%
  mutate(col_split = paste(sample,cluster,treat,sep = "_")) %>%
  split(f = .$col_split) %>%
  lapply(function(x){
    # determine the number of replications of the entry in the metada
    n_row <- x$count
    
    # define the columns
    cluster <- rep(x$cluster,n_row)
    sample <- rep(x$sample,n_row)
    treat <- rep(x$treat,n_row)
    
    # build the fake metadata
    meta <- data.frame(sample,treat,cluster)
    return(meta)
  }) %>%
  bind_rows()

# test stimulation
meta_test <- df %>%
  mutate(col_split = paste(sample,cluster,treat,sep = "_")) %>%
  split(f = .$col_split) %>%
  lapply(function(x){
    # determine the number of replications of the entry in the metada
    n_row <- x$count_test
    
    # define the columns
    cluster <- rep(x$cluster,n_row)
    sample <- rep(x$sample,n_row)
    treat <- rep(x$treat,n_row)
    
    # build the fake metadata
    meta <- data.frame(sample,treat,cluster)
    return(meta)
  }) %>%
  bind_rows()

# check dimensions --------------------------------------------------------
# confirm the numbers from the reference dataset
head(meta_ref)
dim(meta_ref)
sum(df$count)

head(meta_test)
dim(meta_test)
sum(df$count_test)

# run the proportion test diagnosis ---------------------------------------
# Run propeller testing for cell type proportion differences between the groups.
# cluster is the cluster/celltype id
# sample is the reference id of the biological replicates
# group is the grouping id

# ref dataset
properller_out_ref <- propeller(clusters = meta_ref$cluster,
                                sample = meta_ref$sample,
                                group = meta_ref$treat)

properller_out_ref %>%
  rownames_to_column("cell_id") %>%
  write_tsv("../../out/table/propeller_out_ref_simulatedData.tsv")

# test dataset
properller_out_test <- propeller(clusters = meta_test$cluster,
                                 sample = meta_test$sample,
                                 group = meta_test$treat)

properller_out_test %>%
  rownames_to_column("cell_id") %>%
  write_tsv("../../out/table/propeller_out_test_simulatedData.tsv")

# plotting diagnosis ------------------------------------------------------
# reference
df_summary_ref <- meta_ref %>% 
  group_by(cluster,
           sample,
           treat) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(tot = sum(n),
         prop = n/tot) %>%
  ungroup()

# test
df_summary_test <- meta_test %>% 
  group_by(cluster,
           sample,
           treat) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(tot = sum(n),
         prop = n/tot) %>%
  ungroup()

# plot 01
df_summary_ref %>%
  ggplot(aes(x=treat,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~cluster,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/01_propeller_plot01_ref.pdf",width = 10,height = 10)

df_summary_test %>%
  ggplot(aes(x=treat,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~cluster,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/01_propeller_plot01_test.pdf",width = 10,height = 10)

# plot 02
df_summary_ref %>%
  ggplot() +
  geom_boxplot(aes(x=cluster,y=prop,color=treat),outlier.shape = NA) +
  geom_point(aes(x=cluster,y=prop,color=treat),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()+
  ggtitle("summary meta_ref")
ggsave("../../out/plot/01_propeller_plot02_diagnosis_cellid.pdf",width = 8,height = 5)

df_summary_test %>%
  ggplot() +
  geom_boxplot(aes(x=cluster,y=prop,color=treat),outlier.shape = NA) +
  geom_point(aes(x=cluster,y=prop,color=treat),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()+
  ggtitle("summary meta_test")
ggsave("../../out/plot/01_propeller_plot02_test.pdf",width = 8,height = 5)

# plot 03
df_summary_ref %>%
  group_by(treat,cluster) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(treat) %>%
  mutate(tot = sum(n),
         prop = n/tot) %>%
  ggplot() +
  geom_col(aes(x=treat,y=prop,fill=cluster))+
  theme_cowplot()+
  ggtitle("summary meta_ref")

df_summary_test %>%
  group_by(treat,cluster) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(treat) %>%
  mutate(tot = sum(n),
         prop = n/tot) %>%
  ggplot() +
  geom_col(aes(x=treat,y=prop,fill=cluster))+
  theme_cowplot()+
  ggtitle("summary meta_test")
