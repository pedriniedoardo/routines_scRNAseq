# AIM ---------------------------------------------------------------------
# sample analysis using a seurat object

# libraries ---------------------------------------------------------------
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)

# 01 Introduction ---------------------------------------------------------
# Milo is a tool for analysis of complex single cell datasets generated from replicated multi-condition experiments, which detects changes in composition between conditions. While differential abundance (DA) is commonly quantified in discrete cell clusters, Milo uses partially overlapping neighbourhoods of cells on a KNN graph. Starting from a graph that faithfully recapitulates the biology of the cell population, Milo analysis consists of 3 steps:

# Sampling of representative neighbourhoods
# Testing for differential abundance of conditions in all neighbourhoods
# Accounting for multiple hypothesis testing using a weighted FDR procedure that accounts for the overlap of neighbourhoods
# In this vignette we will elaborate on how these steps are implemented in the miloR package.

# 02 Load data ------------------------------------------------------------
data.combined <- readRDS(file = "../../out/object/sobj_total_full_h.rds")
sce <- as.SingleCellExperiment(data.combined)
milo <- Milo(sce)
# add the graph to the object
miloR::graph(milo) <- miloR::graph(buildFromAdjacency(data.combined@graphs$RNA_snn, k=10))
# whit the last command we are adding the graph object calculated in the seurat workflow, in the milo object

# 03 Pre-processing -------------------------------------------------------
# For DA analysis we need to construct an undirected KNN graph of single-cells. Standard single-cell analysis pipelines usually do this from distances in PCA. We normalize and calculate principal components using scater. I also run UMAP for visualization purposes.
# the dataset is already preprocessed
# logcounts(traj_sce) <- log(counts(traj_sce) + 1)
# traj_sce <- runPCA(traj_sce, ncomponents=30)
# traj_sce <- runUMAP(traj_sce)

plotUMAP(milo)

# get some more metadata to plot
colData(milo)

plotUMAP(milo,colour_by="test_cellid")

# 04 Create a Milo object -------------------------------------------------
# For differential abundance analysis on graph neighbourhoods we first construct a Milo object. This extends the SingleCellExperiment class to store information about neighbourhoods on the KNN graph.

# # 4.1From SingleCellExperiment object
# # The Milo constructor takes as input a SingleCellExperiment object.
# traj_milo <- Milo(traj_sce)
# reducedDim(traj_milo, "UMAP") <- reducedDim(traj_sce, "UMAP")
# 
# traj_milo

# # 4.2From AnnData object (.h5ad)
# # We can use the zellkonverter package to make a SingleCellExperiment object from an AnnData object stored as h5ad file.
# library(zellkonverter)
# 
# # Obtaining an example H5AD file.
# example_h5ad <- system.file("extdata", "krumsiek11.h5ad",package = "zellkonverter")
# 
# example_h5ad_sce <- readH5AD(example_h5ad)
# example_h5ad_milo <- Milo(example_h5ad_sce)

# 4.3From Seurat object
# The Seurat package includes a converter to SingleCellExperiment.
# library(Seurat)
# reducedDim(pbmc_small_milo, "PCA", withDimnames=TRUE) <- pbmc_small[['pca']]@cell.embeddings
# reducedDim(pbmc_small_milo, "UMAP", withDimnames=TRUE) <- pbmc_small[['umap']]@cell.embeddings

# 05 Construct KNN graph --------------------------------------------------
# We need to add the KNN graph to the Milo object. This is stored in the graph slot, in igraph format. The miloR package includes functionality to build and store the graph from the PCA dimensions stored in the reducedDim slot.
# skip this if we are usign a seurat object

# traj_milo <- buildGraph(traj_milo, k = 10, d = 30)
# ## Constructing kNN graph with k:10
# # In progress: we are perfecting the functionality to add a precomputed KNN graph (for example constructed with Seurat or scanpy) to the graph slot using the adjacency matrix.
# traj_milo <- buildGraph(traj_milo, k = 10, d = 30)

# 06 1. Defining representative neighbourhoods ----------------------------
# We define the neighbourhood of a cell, the index, as the group of cells connected by an edge in the KNN graph to the index cell. For efficiency, we don’t test for DA in the neighbourhood of every cell, but we sample as indices a subset of representative cells, using a KNN sampling algorithm used by Gut et al. 2015.

# For sampling you need to define a few parameters:

# prop: the proportion of cells to randomly sample to start with (usually 0.1 - 0.2 is sufficient)
# k: the k to use for KNN refinement (we recommend using the same k used for KNN graph building)
# d: the number of reduced dimensions to use for KNN refinement (we recommend using the same d used for KNN graph building)
# refined indicated whether you want to use the sampling refinement algorithm, or just pick cells at random. The default and recommended way to go is to use refinement. The only situation in which you might consider using random instead, is if you have batch corrected your data with a graph based correction algorithm, such as BBKNN, but the results of DA testing will be suboptimal.
milo <- makeNhoods(milo, prop = 0.1, k = 10, d=30, refined = TRUE)
# pbmc_small_milo <- makeNhoods(pbmc_small_milo, prop = 0.1, k = 10, d=30, refined = TRUE)

# Once we have defined neighbourhoods, it’s good to take a look at how big the neighbourhoods are (i.e. how many cells form each neighbourhood). This affects the power of DA testing. We can check this out using the plotNhoodSizeHist function. Empirically, we found it’s best to have a distribution peaking between 50 and 100. Otherwise you might consider rerunning makeNhoods increasing k and/or prop (here the distribution looks ludicrous because it’s a small dataset).
plotNhoodSizeHist(milo)

# 07 Counting cells in neighbourhoods -------------------------------------
# Now we have to count how many cells from each sample are in each neighbourhood. We need to use the cell metadata and specify which column contains the sample information.
milo <- countCells(milo, meta.data = data.frame(colData(milo)), samples="test_donor")
# pbmc_small_milo <- countCells(pbmc_small_milo, meta.data = data.frame(colData(pbmc_small_milo)), samples="groups")

# This adds to the Milo object a n \times m matrix, where n is the number of neighbourhoods and m is the number of experimental samples. Values indicate the number of cells from each sample counted in a neighbourhood. This count matrix will be used for DA testing.
head(nhoodCounts(milo))

# 08 Differential abundance testing ---------------------------------------
# Now we are all set to test for differential abundance in neighbourhoods. We implement this hypothesis testing in a generalized linear model (GLM) framework, specifically using the Negative Binomial GLM implementation in edgeR.

# We first need to think about our experimental design. The design matrix should match samples to a condition of interest. In this case the doxy is the covariate we are going to test for.
traj_design <- data.frame(colData(milo))[,c("test_donor", "test_disease")]

# comparison tested
table(colData(milo)$test_disease)

# comparision in the original dataset
table(data.combined$test_disease)

# potential alternative comparison with more than one level
table(data.combined$test_disease)

traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$test_donor
## Reorder rownames to match columns of nhoodCounts(milo)
traj_design <- traj_design[colnames(nhoodCounts(milo)), , drop=FALSE]

traj_design

# Milo uses an adaptation of the Spatial FDR correction introduced by cydar, which accounts for the overlap between neighbourhoods. Specifically, each hypothesis test P-value is weighted by the reciprocal of the kth nearest neighbour distance. To use this statistic we first need to store the distances between nearest neighbors in the Milo object.
# milo <- calcNhoodDistance(milo, d=30,reduced.dim = "PCA")

# Now we can do the test, explicitly defining our experimental design.
# rownames(traj_design) <- traj_design$sample
da_results <- testNhoods(milo, design = ~test_disease, design.df = traj_design,fdr.weighting = "graph-overlap")

# This calculates a Fold-change and corrected P-value for each neighbourhood, which indicates whether there is significant differential abundance between conditions.
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 

da_results %>%
  arrange(SpatialFDR) %>%
  head()

# 09 Visualize neighbourhoods displaying DA -------------------------------
# To visualize DA results, we build an abstracted graph of neighbourhoods that we can superimpose on the single-cell embedding.
milo <- buildNhoodGraph(milo)

plotUMAP(milo,colour_by="test_cellid") +
  plotNhoodGraphDA(milo, da_results, alpha=1) +
  plot_layout(guides="collect")

# change color of the plot ------------------------------------------------
# saveRDS(milo,"../../out/object/01_miloR.rds")
# saveRDS(da_results,"../../out/object/01_da_results.rds")
# milo <- readRDS("out/object/miloR.rds")
# da_results <- readRDS("out/object/da_results.rds")

#
test <- plotNhoodGraphDA(milo, da_results, alpha=1)
test+scale_fill_viridis_c(option = "turbo")

test+scale_fill_gradient2(
  low = "blue",
  mid = "white",
  high = "red",
  midpoint = 0,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill"
)

test$data %>% 
  ggplot(aes(x=x,y=y))+
  geom_point(aes(fill=colour_by,size = size),shape=21,alpha=0.5)+scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  )+theme_bw()
ggsave("../../out/plot/01_miloR_test_update2.pdf",width = 7,height = 6)

# USING A SEURAT INTEGRATED OBJECT
# the general workflow remains the same, the main difference is the fact that at the beginnig, when we are generating the milo object, we should use the following:
sce <- as.SingleCellExperiment(data.combined,assay = "RNA")
milo <- Milo(sce)

# add the graph to the object
miloR::graph(milo) <- miloR::graph(buildFromAdjacency(data.combined@graphs$integrated_snn, k=10))

# MORE THAN ONE COMPARISON
# in some situations is important to be able to make more than one comparison from the dataset. There is a specific sample case for this situation:

# libraries ---------------------------------------------------------------
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(MouseThymusAgeing)
library(scuttle)

# Introduction ------------------------------------------------------------
# We have seen how Milo uses graph neighbourhoods to model cell state abundance differences in an experiment, when comparing 2 groups. However, we are often interested in testing between 2 specific groups in our analysis when our experiment has collected data from $\gt$  2 groups. We can focus our analysis to a 2 group comparison and still make use of all of the data for things like dispersion estimation, by using _contrasts_. For an in-depth use of contrasts we recommend users refer to the `limma`  and `edgeR` Biostars and Bioconductor community forum threads on the subject. Here I will give an overview of how to use contrasts in the context of a Milo analysis.

# Load data ---------------------------------------------------------------
# We will use the `MouseThymusAgeing` data package as there are multiple groups that we can compare.
thy.sce <- MouseSMARTseqData() # this function downloads the full SCE object
thy.sce <- logNormCounts(thy.sce)
thy.sce

# Define cell neighbourhoods ----------------------------------------------
thy.sce <- runUMAP(thy.sce) # add a UMAP for plotting results later
thy.milo <- Milo(thy.sce) # from the SCE object
reducedDim(thy.milo, "UMAP") <- reducedDim(thy.sce, "UMAP")
plotUMAP(thy.milo, colour_by="SubType") + plotUMAP(thy.milo, colour_by="Age")

# These UMAPs shows how the different thymic epithelial cell subtypes and cells from different aged mice are distributed across our single-cell data set. Next we build the KNN graph and define neighbourhoods to quantify cell abundance across our experimental samples.

# build KNN graph ---------------------------------------------------------
thy.milo <- buildGraph(thy.milo, k = 10, d = 30)
thy.milo <- makeNhoods(thy.milo, prop = 0.1, k = 10, d=30, refined = TRUE, refinement_scheme="graph") # make nhoods using graph-only as this is faster
colData(thy.milo)$Sample <- paste(colData(thy.milo)$SortDay, colData(thy.milo)$Age, sep="_")
thy.milo <- countCells(thy.milo, meta.data = data.frame(colData(thy.milo)), samples="Sample") # make the nhood X sample counts matrix

# Differential abundance testing with contrasts ---------------------------
# Now we have the pieces in place for DA testing to demonstrate how to use contrasts. We will use these contrasts to explicitly define which groups will be compared to each other.
thy.design <- data.frame(colData(thy.milo))[,c("Sample", "SortDay", "Age")]
thy.design <- distinct(thy.design)
rownames(thy.design) <- thy.design$Sample
## Reorder rownames to match columns of nhoodCounts(milo)
thy.design <- thy.design[colnames(nhoodCounts(thy.milo)), , drop=FALSE]

thy.design

# see the number of different levels in the model
table(thy.design$Age)

# To demonstrate the use of contrasts we will fit the whole model to the whole data set, but we will compare sequential pairs of time points. I'll start with week 1 vs.week 4 to illustrate the syntax.
rownames(thy.design) <- thy.design$Sample
contrast.1 <- c("Age1wk - Age4wk") # the syntax is <VariableName><ConditionLevel> - <VariableName><ControlLevel>
# we need to use the ~ 0 + Variable expression here so that we have all of the levels of our variable as separate columns in our model matrix
da_results <- testNhoods(thy.milo, design = ~ 0 + Age, design.df = thy.design, model.contrasts = contrast.1, fdr.weighting="graph-overlap")
head(da_results)

table(da_results$SpatialFDR < 0.1)

# This calculates a Fold-change and corrected P-value for each neighbourhood, which indicates whether there is significant differential abundance between conditions for `r sum(da_results$SpatialFDR < 0.1)` neighbourhoods.

# You will notice that the syntax for the contrasts is quite specific. It starts with the name of the column variable that contains the different group levels; in this case it is the `Age` variable. We then define the comparison levels as `level1 - level2`. To understand this syntax we need to consider what we are concretely comparing. In this case we are asking what is the ratio of the average cell count at week1 compared to the average cell count at week 4, where the averaging is across the replicates. The  reason we express this as a difference rather than a ratio is because we are dealing with the _log_ fold change.

# We can also pass multiple comparisons at the same time, for instance if we wished to compare each sequential pair of time points. This will give us a better intuition behind how to use contrasts to compare multiple groups.

contrast.all <- c("Age1wk - Age4wk", "Age4wk - Age16wk", "Age16wk - Age32wk", "Age32wk - Age52wk")
# this is the edgeR code called by `testNhoods`
model <- model.matrix(~ 0 + Age, data=thy.design)
mod.constrast <- makeContrasts(contrasts=contrast.all, levels=model)
mod.constrast

list_da <- lapply(contrast.all,function(x){
  testNhoods(thy.milo, design = ~ 0 + Age, design.df = thy.design, model.contrasts = x, fdr.weighting="graph-overlap")
}) %>% 
  setNames(contrast.all)

lapply(list_da,function(x){
  head(x)
})

# This shows the contrast matrix. If we want to test each of these comparisons then we need to pass them sequentially to `testNhoods`, then apply an additional  multiple testing correction to the spatial FDR values.

# Contrasts are not limited to these simple pair-wise comparisons, we can also group levels together for comparisons. For instance, imagine we want to know what the effect of the cell counts in the week 1 mice is _compared to all other time points_.
model <- model.matrix(~ 0 + Age, data=thy.design)
ave.contrast <- c("Age1wk - (Age4wk + Age16wk + Age32wk + Age52wk)/4")
mod.constrast <- makeContrasts(contrasts=ave.contrast, levels=model)
mod.constrast

# In this contrasts matrix we can see that we have taken the average effect over the other time points. Now running this using `testNhoods`
da_results <- testNhoods(thy.milo, design = ~ 0 + Age, design.df = thy.design, model.contrasts = ave.contrast, fdr.weighting="graph-overlap")
table(da_results$SpatialFDR < 0.1)

# In this comparison there are `r sum(da_results$SpatialFDR < 0.1)` DA nhoods - which we can visualise on a superimposed single-cell UMAP.
thy.milo <- buildNhoodGraph(thy.milo)
plotUMAP(thy.milo, colour_by="SubType") + plotNhoodGraphDA(thy.milo, da_results, alpha=0.1) +
  plot_layout(guides="auto" )

# In these side-by-side UMAPs we can see that there is an enrichment of the Perinatal cTEC and Proliferating TEC populations in the 1 week old compared tothe other time points.

# For a more extensive description of the uses of contrasts please take a look at the edgeR documentation \Biocpkg{edgeR}.
