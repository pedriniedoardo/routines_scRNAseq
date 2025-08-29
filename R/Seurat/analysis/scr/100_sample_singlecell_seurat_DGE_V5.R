# AIM ---------------------------------------------------------------------
# sample routine for DE analysis at the single cell level with multiple levels

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(limma)
library(ggrepel)
library(presto)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
data.combined <- readRDS("../out/object/100_ifnb_DonorStim.rds")

# check the object version
class(data.combined@assays$RNA)

# compare the clusters per treatment --------------------------------------
# make suer the correct defalult dataset il loaded should be RNA
DefaultAssay(data.combined)

# check that in the RNA slot the data object is indeed loaded with normalized values
data.combined@assays$RNA$data[1:10,1:10]

# define the grouping variables for the comparison of stim
head(data.combined@meta.data)

# define the new grouping
table(data.combined@meta.data$seurat_annotations,data.combined@meta.data$stim)

data.combined$condition.annotation <- paste(data.combined$stim,data.combined$seurat_annotations, sep = "_")
head(data.combined@meta.data)

# update the idents of the object
Idents(data.combined) <- "condition.annotation"

# avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
# run the differential expression over the WT
# clusters_id <- as.character(sort(unique(data.combined$annotation_confident))) %>%
#   str_subset(pattern = "P09",negate = T)
clusters_id <- as.character(sort(unique(data.combined$seurat_annotations)))

# NOLD vs control all
# run presto to make the comparison faster
list_STIM_vs_CTRL <- lapply(clusters_id,function(x){
  id_1 <- paste0("STIM_",x)
  id_2 <- paste0("CTRL_",x)
  # latest version of seurat might have implemented the presto methods by default, therefore the regular FindMarkers function might work well
  # response <- FindMarkers(data.combined, ident.1 = id_1, ident.2 = id_2, verbose = T,logfc.threshold = 0)
  response <- RunPresto(data.combined, ident.1 = id_1, ident.2 = id_2, verbose = T,logfc.threshold = 0)
  response %>%
    rownames_to_column("gene") %>%
    mutate(id_1 = id_1,
           id_2 = id_2) %>%
    mutate(cluster = x)
})

list_STIM_vs_CTRL %>%
  setNames(paste0(clusters_id,"_STIM_vs_CTRL")) %>%
  bind_rows(.id = "annotation") %>%
  write_tsv("../out/table/100_response_STIM_vs_CTRL_seurat_sc_V5.tsv")

# plot the volcanos per cluster -------------------------------------------
folder <- "../out/table/"
file <- dir(folder) %>%
  str_subset(pattern = "100_response_") %>% 
  str_subset(pattern = "STIM_vs_CTRL") %>%
  str_subset(pattern = "V5")

df_res <- lapply(file, function(x){
  test_plot <- read_tsv(paste0(folder,x))
}) %>%
  bind_rows()

# show the distribution of the pvalues
df_res %>%
  mutate(comparison = str_extract(annotation,pattern = "STIM_vs_CTRL")) %>%
  ggplot(aes(x = p_val)) + geom_histogram()+theme_bw()+facet_grid(comparison~cluster)+theme(strip.background = element_blank())
ggsave("../out/plot/100_dist_p_value_STIM_vs_CTRL_seurat_sc_V5.pdf",width = 20,height = 4)

# render all of them as a volcano plot
test_significant <- df_res %>%
  mutate(comparison = str_extract(annotation,pattern = "STIM_vs_CTRL")) %>%
  mutate(threshold = case_when(abs(avg_log2FC) > 1 & p_val_adj<0.05~1,
                               T~0)) %>%
  filter(threshold == 1)

# library(ggrepel)
df_res %>%
  mutate(comparison = str_extract(annotation,pattern = "STIM_vs_CTRL")) %>%
  # filter(symbol %in% setdiff(GOI_SLC,GOI)) %>%
  ggplot(aes(x = avg_log2FC,y = -log(p_val_adj))) +
  geom_point(alpha = 0.01) +
  geom_point(data = test_significant,aes(x = avg_log2FC,y = -log(p_val_adj)),col="red",alpha = 0.5) +
  geom_vline(xintercept = c(-1,1),col="gray",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="gray",linetype="dashed")+
  geom_text_repel(data = test_significant,aes(x = avg_log2FC,y = -log(p_val_adj),label = gene)) +
  facet_grid(comparison~cluster) +
  theme_bw()+theme(strip.background = element_blank())
ggsave("../out/plot/100_volcano_STIM_vs_CTRL_seurat_sc_V5.pdf",width = 20,height = 5)
