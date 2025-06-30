# AIM ---------------------------------------------------------------------
# load the metadata from the initial processing and produce the QC plots

# libraries ---------------------------------------------------------------
library(tidyverse)

# read in the data --------------------------------------------------------
file <- dir("../out/table",full.names = T) %>%
  str_subset(pattern = "_meta_preQC.tsv")

sample_name <- str_remove(file,pattern = "../out/table/01_") %>%
  str_remove(pattern = "_meta_preQC.tsv")

meta_total <- pmap(list(file,sample_name),function(x,x_name){
  read_tsv(x) %>%
    mutate(dataset = x_name)
}) %>%
  bind_rows()

# define some plotting parameters
# define number of rows in the panel
test <- length(sample_name)
# test <- 10
panel_row <- test %>% sqrt() %>% round(digits = 0)
panel_col <- (test / panel_row) %>% ceiling()

# plot QC -----------------------------------------------------------------
# save the total metadata before QC
meta_total %>%
  write_tsv("../out/table/meta_total_beforeQC_V5.tsv")

# fixed threshold scatter nFeature vs percent.mt
meta_total %>%
  ggplot(aes(y = percent.mt,x = nFeature_RNA)) + geom_point(alpha=0.3) +
  facet_wrap(~orig.ident) +
  theme_bw() +
  theme(strip.background = element_blank())
# save the plot
ggsave("../out/plot/01_fixed_scatter_feature_mito_V5.pdf",
       width = 4 * panel_col,
       height = 4 * panel_row)

# adaptive threshold scatter nFeature vs percent.mt single value
meta_total %>%
  ggplot(aes(y = percent.mt,x = nFeature_RNA,col=discard_single)) + geom_point(alpha=0.3) +
  facet_wrap(~orig.ident) +
  theme_bw() +
  theme(strip.background = element_blank())
# save the plot
ggsave("../out/plot/01_fixed_scatter_feature_mito_automaticSingle_V5.pdf",
       width = 4 * panel_col,
       height = 4 * panel_row)

# adaptive threshold scatter nFeature vs percent.mt multi value
meta_total %>%
  ggplot(aes(y = percent.mt,x = nFeature_RNA,col=discard_multi)) + geom_point(alpha=0.3) +
  facet_wrap(~orig.ident) +
  theme_bw() +
  theme(strip.background = element_blank())
# save the plot
ggsave("../out/plot/01_fixed_scatter_feature_mito_automaticMulti_V5.pdf",
       width = 4 * panel_col,
       height = 4 * panel_row)

#
meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  ggplot(aes(x=orig.ident,y=value)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.01) +
  facet_wrap(~var,scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  theme(strip.background = element_blank())

# save the plot
ggsave("../out/plot/01_fixed_boxplot_reads_V5.pdf",
       width = 2 * test + 2,
       height = 8)

#
meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "percent.mt") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
  facet_wrap(orig.ident~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(0,1,2,5,10,20,40,60,80,100)) +
  geom_vline(xintercept = c(10),col="red",linetype="dashed") +
  annotate("rect", xmin=0, xmax=10, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())

# save the plot
ggsave("../out/plot/01_fixed_histo_mito_V5.pdf",
       width = 4 * panel_col,
       height = 4 * panel_row)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nFeature_RNA") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(orig.ident~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(200,500,1000,2000,4000,6000,8000,10000,20000)) +
  geom_vline(xintercept = c(1000,6000),col="red",linetype="dashed") +
  annotate("rect", xmin=1000, xmax=6000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# save the plot
ggsave("../out/plot/01_fixed_histo_features_V5.pdf",
       width = 4 * panel_col,
       height = 4 * panel_row)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nCount_RNA") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(orig.ident~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(100,500,1000,2000,4000,8000,20000,40000,80000,200000)) +
  # geom_vline(xintercept = c(500,5000),col="red",linetype="dashed") +
  # annotate("rect", xmin=500, xmax=5000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# save the plot
ggsave("../out/plot/01_fixed_histo_counts_V5.pdf",
       width = 4 * panel_col,
       height = 4 * panel_row)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.ribo,percent.globin)) %>%
  filter(var == "percent.ribo") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(orig.ident~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(0,1,2,5,10,20,40,60,80,100)) +
  # geom_vline(xintercept = c(500,5000),col="red",linetype="dashed") +
  # annotate("rect", xmin=500, xmax=5000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# save the plot
ggsave("../out/plot/01_fixed_histo_ribo_V5.pdf",
       width = 4 * panel_col,
       height = 4 * panel_row)

meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.ribo,percent.globin)) %>%
  filter(var == "percent.globin") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(orig.ident~var,scales = "free") +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(0,1,2,5,10,20,40,60,80,100)) +
  # geom_vline(xintercept = c(500,5000),col="red",linetype="dashed") +
  # annotate("rect", xmin=500, xmax=5000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# save the plot
ggsave("../out/plot/01_fixed_histo_globin_V5.pdf",
       width = 4 * panel_col,
       height = 4 * panel_row)
