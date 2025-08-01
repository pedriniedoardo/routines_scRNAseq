# AIM ---------------------------------------------------------------------
# they idea is to try to gauge what kind of information is acquired from including or not the introns in the analysis.
# the first approach is to load the matrices generated by w introns or wo introns. make sure we are using the same cells and  aggregare the overall reads per gene. the difference in reads is accounted by intronic reads only. assign the proportion of intronic reads per gene.

# libraries ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(DoubletFinder)
library(UpSetR)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# define the filtering parameters -----------------------------------------
featureLow_thr <- 1000
featureHigh_thr <- 6000
mito_thr <- 15
label <- "01000_06000_15_V5"

# read in the data --------------------------------------------------------
# location of all the raw matrices
id_sample <- dir("../data/test_introns/cellranger7_out/")

data_w_introns <- Read10X(data.dir = paste0("../data/test_introns/cellranger7_out/sample_untreated_Wintron/filtered_feature_bc_matrix/"))
data_wo_introns <- Read10X(data.dir = paste0("../data/test_introns/cellranger7_out/sample_untreated_WOintron/filtered_feature_bc_matrix/"))

# wrangling ---------------------------------------------------------------
# make sure to compare the same cells in both dataset
# how many cells are in common
test_intron <- list(barcodes_w_intron = colnames(data_w_introns),
     barcodes_wo_intron = colnames(data_wo_introns))

UpSetR::upset(fromList(test_intron), nsets = 2, sets.bar.color = "#56B4E9", order.by = "freq")

# focus only on the common barcodes
common_cells <- intersect(colnames(data_w_introns),colnames(data_wo_introns))

# make sure the genes are matched. if so the difference can be done directly
sum(!rownames(data_w_introns) == rownames(data_wo_introns))

# compunds the reads per gene with both version
sobj_common_w_intron <- CreateSeuratObject(counts = data_w_introns[,common_cells],
                                           min.cells = 0,
                                           min.features = 0)

sobj_common_wo_intron <- CreateSeuratObject(counts = data_wo_introns[,common_cells],
                                           min.cells = 0,
                                           min.features = 0)


# maks sure the gene and the colnames are the same
sum(!rownames(sobj_common_w_intron) == rownames(sobj_common_wo_intron))
sum(!colnames(sobj_common_w_intron) == colnames(sobj_common_wo_intron))

# aggregate the reads per gene in both modalities
df_count_w_intron <- AggregateExpression(sobj_common_w_intron)$RNA %>%
  data.frame() %>%
  setNames(c("count_w_intron")) %>%
  rownames_to_column("gene")

df_count_wo_intron <- AggregateExpression(sobj_common_wo_intron)$RNA %>%
  data.frame() %>%
  setNames(c("count_wo_intron")) %>%
  rownames_to_column("gene")

# pull some metadata from the genes

# join the two tables
df_summary <- left_join(df_count_w_intron,df_count_wo_intron,by = "gene") %>%
  # add a pseudocount of 1 to all the genes
  mutate(count_w_intron_adj = count_w_intron + 1,
         count_wo_intron_adj = count_wo_intron + 1) %>%
  mutate(diff = count_w_intron - count_wo_intron)

# scatter plot w vs wo
p1 <- df_summary %>%
  ggplot(aes(x = count_wo_intron_adj, y = count_w_intron_adj)) +
  geom_point(shape = 1,
             alpha=0.3) +
  geom_abline(intercept = 0,slope = 1,col="red",linetype="dashed") +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10()

# plot the delta
df_summary %>%
  filter(diff < 0)
