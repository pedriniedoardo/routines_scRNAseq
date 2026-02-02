# # AIM ---------------------------------------------------------------------
# # test using the dataset shared by Sofia
# 
# # libraries ---------------------------------------------------------------
# library(miloR)
# library(SingleCellExperiment)
# library(scater)
# library(scran)
# library(dplyr)
# library(patchwork)
# library(Seurat)
# 
# # 01 Introduction ---------------------------------------------------------
# # Milo is a tool for analysis of complex single cell datasets generated from replicated multi-condition experiments, which detects changes in composition between conditions. While differential abundance (DA) is commonly quantified in discrete cell clusters, Milo uses partially overlapping neighbourhoods of cells on a KNN graph. Starting from a graph that faithfully recapitulates the biology of the cell population, Milo analysis consists of 3 steps:
# 
# # Sampling of representative neighbourhoods
# # Testing for differential abundance of conditions in all neighbourhoods
# # Accounting for multiple hypothesis testing using a weighted FDR procedure that accounts for the overlap of neighbourhoods
# # In this vignette we will elaborate on how these steps are implemented in the milo3R package.
# 
# # 02 Load data ------------------------------------------------------------
# ref <- readRDS(file = "../../data/harmonyHO.rds")
# # c("Ctr 5 day","FOP 5 day","Ctr 7 day","FOP 7 day")
# # test <- subset(ref,subset = condition %in% c("Ctr 5 day","FOP 5 day","Ctr 7 day","FOP 7 day"))
# 
# DimPlot(ref)
# class(ref@assays$RNA)
# 
# sce3 <- as.SingleCellExperiment(ref)
# milo3 <- Milo(sce3)
# 
# DefaultAssay(ref)
# 
# str(ref@graphs$RNA_snn)
# str(ref@graphs$RNA_nn)
# 
# # buildFromAdjacency(test@graphs$RNA_snn, k=10)
# # add the graph to the object
# miloR::graph(milo3) <- miloR::graph(buildFromAdjacency(ref@graphs$RNA_snn, k=10,is.binary = F))
# # whit the last command we are adding the graph object calculated in the seurat workflow, in the milo3 object
# 
# # 03 Pre-processing -------------------------------------------------------
# # For DA analysis we need to construct an undirected KNN graph of single-cells. Standard single-cell analysis pipelines usually do this from distances in PCA. We normalize and calculate principal components using scater. I also run UMAP for visualization purposes.
# # the dataset is already preprocessed
# # logcounts(traj_sce2) <- log(counts(traj_sce2) + 1)
# # traj_sce2 <- runPCA(traj_sce2, ncomponents=30)
# # traj_sce2 <- runUMAP(traj_sce2)
# 
# plotUMAP(milo3)
# 
# # get some more metadata to plot
# colData(milo3)
# 
# plotUMAP(milo3,colour_by="ident")
# 
# # 04 Create a Milo object -------------------------------------------------
# # For differential abundance analysis on graph neighbourhoods we first construct a Milo object. This extends the SingleCellExperiment class to store information about neighbourhoods on the KNN graph.
# 
# # # 4.1From SingleCellExperiment object
# # # The Milo constructor takes as input a SingleCellExperiment object.
# # traj_milo3 <- Milo(traj_sce2)
# # reducedDim(traj_milo3, "UMAP") <- reducedDim(traj_sce2, "UMAP")
# # 
# # traj_milo3
# 
# # # 4.2From AnnData object (.h5ad)
# # # We can use the zellkonverter package to make a SingleCellExperiment object from an AnnData object stored as h5ad file.
# # library(zellkonverter)
# # 
# # # Obtaining an example H5AD file.
# # example_h5ad <- system.file("extdata", "krumsiek11.h5ad",package = "zellkonverter")
# # 
# # example_h5ad_sce2 <- readH5AD(example_h5ad)
# # example_h5ad_milo3 <- Milo(example_h5ad_sce2)
# 
# # 4.3From Seurat object
# # The Seurat package includes a converter to SingleCellExperiment.
# # library(Seurat)
# # reducedDim(pbmc_small_milo3, "PCA", withDimnames=TRUE) <- pbmc_small[['pca']]@cell.embeddings
# # reducedDim(pbmc_small_milo3, "UMAP", withDimnames=TRUE) <- pbmc_small[['umap']]@cell.embeddings
# 
# # 05 Construct KNN graph --------------------------------------------------
# # We need to add the KNN graph to the Milo object. This is stored in the graph slot, in igraph format. The milo3R package includes functionality to build and store the graph from the PCA dimensions stored in the reducedDim slot.
# # skip this if we are usign a seurat object
# 
# # traj_milo3 <- buildGraph(traj_milo3, k = 10, d = 30)
# # ## Constructing kNN graph with k:10
# # # In progress: we are perfecting the functionality to add a precomputed KNN graph (for example constructed with Seurat or scanpy) to the graph slot using the adjacency matrix.
# # traj_milo3 <- buildGraph(traj_milo3, k = 10, d = 30)
# 
# # 06 1. Defining representative neighbourhoods ----------------------------
# # We define the neighbourhood of a cell, the index, as the group of cells connected by an edge in the KNN graph to the index cell. For efficiency, we don’t test for DA in the neighbourhood of every cell, but we sample as indices a subset of representative cells, using a KNN sampling algorithm used by Gut et al. 2015.
# 
# # For sampling you need to define a few parameters:
# 
# # prop: the proportion of cells to randomly sample to start with (usually 0.1 - 0.2 is sufficient)
# # k: the k to use for KNN refinement (we recommend using the same k used for KNN graph building)
# # d: the number of reduced dimensions to use for KNN refinement (we recommend using the same d used for KNN graph building)
# # refined indicated whether you want to use the sampling refinement algorithm, or just pick cells at random. The default and recommended way to go is to use refinement. The only situation in which you might consider using random instead, is if you have batch corrected your data with a graph based correction algorithm, such as BBKNN, but the results of DA testing will be suboptimal.
# milo3 <- makeNhoods(milo3, prop = 0.1, k = 10, d=30, refined = TRUE)
# # pbmc_small_milo3 <- makeNhoods(pbmc_small_milo3, prop = 0.1, k = 10, d=30, refined = TRUE)
# 
# # Once we have defined neighbourhoods, it’s good to take a look at how big the neighbourhoods are (i.e. how many cells form each neighbourhood). This affects the power of DA testing. We can check this out using the plotNhoodSizeHist function. Empirically, we found it’s best to have a distribution peaking between 50 and 100. Otherwise you might consider rerunning makeNhoods increasing k and/or prop (here the distribution looks ludicrous because it’s a small dataset).
# plotNhoodSizeHist(milo3)
# 
# # 07 Counting cells in neighbourhoods -------------------------------------
# # Now we have to count how many cells from each sample are in each neighbourhood. We need to use the cell metadata and specify which column contains the sample information.
# milo3 <- countCells(milo3, meta.data = data.frame(colData(milo3)), samples="test_donor")
# # pbmc_small_milo3 <- countCells(pbmc_small_milo3, meta.data = data.frame(colData(pbmc_small_milo3)), samples="groups")
# 
# # This adds to the Milo object a n \times m matrix, where n is the number of neighbourhoods and m is the number of experimental samples. Values indicate the number of cells from each sample counted in a neighbourhood. This count matrix will be used for DA testing.
# head(nhoodCounts(milo3))
# 
# # 08 Differential abundance testing ---------------------------------------
# # Now we are all set to test for differential abundance in neighbourhoods. We implement this hypothesis testing in a generalized linear model (GLM) framework, specifically using the Negative Binomial GLM implementation in edgeR.
# 
# # We first need to think about our experimental design. The design matrix should match samples to a condition of interest. In this case the doxy is the covariate we are going to test for.
# traj_design2 <- data.frame(colData(milo3))[,c("test_donor", "test_disease")]
# 
# # comparison tested
# table(colData(milo3)$test_disease)
# 
# # comparision in the original dataset
# table(test$test_disease)
# 
# # potential alternative comparison with more than one level
# table(test$test_disease)
# 
# traj_design2 <- distinct(traj_design2)
# rownames(traj_design2) <- traj_design2$test_donor
# ## Reorder rownames to match columns of nhoodCounts(milo3)
# traj_design2 <- traj_design2[colnames(nhoodCounts(milo3)), , drop=FALSE]
# 
# traj_design2
# 
# # Milo uses an adaptation of the Spatial FDR correction introduced by cydar, which accounts for the overlap between neighbourhoods. Specifically, each hypothesis test P-value is weighted by the reciprocal of the kth nearest neighbour distance. To use this statistic we first need to store the distances between nearest neighbors in the Milo object.
# # milo3 <- calcNhoodDistance(milo3, d=30,reduced.dim = "PCA")
# 
# # Now we can do the test, explicitly defining our experimental design.
# # rownames(traj_design2) <- traj_design2$sample
# da_results2 <- testNhoods(milo3, design = ~test_disease, design.df = traj_design2,fdr.weighting = "graph-overlap")
# 
# # This calculates a Fold-change and corrected P-value for each neighbourhood, which indicates whether there is significant differential abundance between conditions.
# da_results2 %>%
#   arrange(- SpatialFDR) %>%
#   head() 
# 
# da_results2 %>%
#   arrange(SpatialFDR) %>%
#   head()
# 
# # 09 Visualize neighbourhoods displaying DA -------------------------------
# # To visualize DA results, we build an abstracted graph of neighbourhoods that we can superimpose on the single-cell embedding.
# milo3 <- buildNhoodGraph(milo3)
# 
# plotUMAP(milo3,colour_by="test_cellid") +
#   plotNhoodGraphDA(milo3, da_results2, alpha=1) +
#   plot_layout(guides="collect")
# 
# # change color of the plot ------------------------------------------------
# #
# test2 <- plotNhoodGraphDA(milo3, da_results2, alpha=1)
# test2+scale_fill_viridis_c(option = "turbo")
# 
# test2+scale_fill_gradient2(
#   low = "blue",
#   mid = "white",
#   high = "red",
#   midpoint = 0,
#   space = "Lab",
#   na.value = "grey50",
#   guide = "colourbar",
#   aesthetics = "fill"
# )
# 
# test2$data %>% 
#   ggplot(aes(x=x,y=y))+
#   geom_point(aes(fill=colour_by,size = size),shape=21,alpha=0.5)+scale_fill_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     midpoint = 0,
#     space = "Lab",
#     na.value = "grey50",
#     guide = "colourbar",
#     aesthetics = "fill"
#   )+theme_bw()
# ggsave("out/image/milo3R_test_update2.pdf",width = 7,height = 6)
# 
# # Inspecting DA testing results -------------------------------------------
# # We can start inspecting the results of our DA analysis from a couple of standard diagnostic plots. We first inspect the distribution of uncorrected P values, to verify that the test was balanced.
# ggplot(da_results2, aes(PValue)) + geom_histogram(bins=50)+theme_bw()
# 
# ggplot(da_results2, aes(logFC, -log10(SpatialFDR))) + 
#   geom_point(shape = 1,alpha=0.5) +
#   ## Mark significance threshold (5% FDR)
#   geom_hline(yintercept = -log10(0.05),linetype="dashed",col="gray") +
#   theme_bw()+
#   theme(strip.background = element_blank())
# 
# # To visualize DA results relating them to the embedding of single cells, we can build an abstracted graph of neighbourhoods that we can superimpose on the single-cell embedding. Here each node represents a neighbourhood, while edges indicate how many cells two neighbourhoods have in common. Here the layout of nodes is determined by the position of the index cell in the UMAP embedding of all single-cells. The neighbourhoods displaying significant DA are colored by their log-Fold Change.
# milo3 <- buildNhoodGraph(milo3)
# 
# ## Plot single-cell UMAP
# umap_pl2 <- plotReducedDim(milo3, dimred = "UMAP", colour_by="test_disease", text_by = "test_cellid", 
#                            text_size = 3, point_size=0.5) +
#   guides(fill="none")
# 
# ## Plot neighbourhood graph
# nh_graph_pl2 <- plotNhoodGraphDA(milo3, da_results2, layout="UMAP",alpha=0.98) 
# umap_pl2 + nh_graph_pl2 + plot_layout(guides="collect")
# 
# # We might also be interested in visualizing whether DA is particularly evident in certain cell types. To do this, we assign a cell type label to each neighbourhood by finding the most abundant cell type within cells in each neighbourhood. We can label neighbourhoods in the results data.frame using the function annotateNhoods. This also saves the fraction of cells harbouring the label.
# da_results2 <- annotateNhoods(milo3, da_results2, coldata_col = "test_cellid")
# head(da_results2)
# 
# # Then we visualize the test results with a volcano plot (remember that each point here represents a neighbourhood, not a cell).
# ggplot(da_results2, aes(logFC, -log10(SpatialFDR))) + 
#   geom_point(shape = 1,alpha=0.5) +
#   ## Mark significance threshold (5% FDR)
#   geom_hline(yintercept = -log10(0.05),linetype="dashed",col="gray") +
#   theme_bw()+
#   theme(strip.background = element_blank())+
#   facet_wrap(~test_cellid)
# 
# # While neighbourhoods tend to be homogeneous, we can define a threshold for celltype_fraction to exclude neighbourhoods that are a mix of cell types.
# ggplot(da_results2, aes(test_cellid_fraction)) + geom_histogram(bins=50)
# 
# da_results2$test_cellid <- ifelse(da_results2$test_cellid_fraction < 0.7, "Mixed", da_results2$test_cellid)
# # Now we can visualize the distribution of DA Fold Changes in different cell types
# test2 <- plotDAbeeswarm(da_results2, group.by = "test_cellid",alpha = 1)
# 
# plotDAbeeswarm
# 
# test2$data %>%
#   group_by(group_by) %>%
#   summarise(med = median(pos_x)) %>%
#   arrange(med) %>%
#   pull(group_by)
# 
# test2$data %>%
#   mutate(color = case_when(FDR<0.05~"sig",
#                            T~"n.s.")) %>%
#   arrange(color) %>%
#   ggplot(aes(x=logFC,y=group_by,color = color)) +
#   ggbeeswarm::geom_quasirandom(alpha=0.5)+
#   theme_bw()+
#   scale_color_manual(values = c("gray","red"))+
#   scale_x_continuous(limits = c(-2.5,2.5))
# 
# # This is already quite informative: we can see that certain early development cell types, such as epiblast and primitive streak, are enriched in the earliest time stage, while others are enriched later in development, such as ectoderm cells. Interestingly, we also see plenty of DA neighbourhood with a mixed label. This could indicate that transitional states show changes in abundance in time.
