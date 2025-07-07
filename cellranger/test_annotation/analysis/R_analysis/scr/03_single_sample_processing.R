# AIM ---------------------------------------------------------------------
# single standard sample preprocessing

# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)
library(homologene)
library(limma)
library(Polychrome)
library(pals)
library(presto)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1000 * 1024^2)

# read in the data --------------------------------------------------------
# read in the filtered dataset
list_datasc <- readRDS("../../out/R_analysis/object/02_list_datasc_SoupX_afterQC_1000_7000_10_V5_singlet.rds")

# preprocessing -----------------------------------------------------------
# x <- list_datasc$test_neuron_auto
list_datasc_pro <- lapply(list_datasc,function(x){
  # track the processing
  print(x)
  DefaultAssay(x) <- "RNA"
  x@assays$RNA$counts[1:10,1:10]
  
  # normalize the objet
  x <- x %>%
    NormalizeData(verbose = T)
  
  # convert them to mouse gene
  # library(homologene)
  # s.genes <- homologene(cc.genes$s.genes, inTax = 9606, outTax = 10090) %>%
  #   pull("10090")
  # g2m.genes <- homologene(cc.genes$g2m.genes, inTax = 9606, outTax = 10090) %>%
  #   pull("10090")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes)
  
  x %>%
    # skip the normalizatio that has been already performed at the beginning
    # Seurat::NormalizeData(verbose = T) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>%
    # I can scale the missing features afterwards now focus on the highly variable one for speed purposes
    # ScaleData(vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"), verbose = T) %>% 
    ScaleData(vars.to.regress = c("percent.mt","nCount_RNA"), verbose = T) %>% 
    # run this if you want to scale all the variables
    # ScaleData(vars.to.regress = c("percent.mt.harmony","nCount_RNA.harmony","S.Score.harmony","G2M.Score.harmony"), verbose = T,features = all.genes) %>% 
    RunPCA(npcs = 30, verbose = T) %>% 
    RunUMAP(reduction = "pca", dims = 1:30,return.model = TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
    identity()
})
# save the object post processing
saveRDS(list_datasc_pro,"../../out/R_analysis/object/03_list_datasc_fixed_pro.rds")

# plotting for the report -------------------------------------------------
# list_datasc_pro <- readRDS("../../out/R_analysis/object/03_list_datasc_fixed_pro.rds")

# loop the plotting for each object
# x <- list_datasc_pro$SacsKO5mo_MG
lapply(list_datasc_pro,function(x){
  id_sample <- unique(x@meta.data$orig.ident)
  
  # save the array of resolutions -------------------------------------------
  id_resolution <- str_subset(colnames(x@meta.data),pattern = "RNA_snn_res") %>%
    sort()
  
  list_plot <- lapply(id_resolution,function(res){
    plot <- DimPlot(x,
                    reduction = "umap",
                    group.by = res,
                    label = T,
                    raster = F)
    return(plot)
  })
  
  wrap_plots(list_plot)
  ggsave(paste0("../../out/R_analysis/plot/UMAPCluster_resolutions_",id_sample,".pdf"),width = 25,height = 15)
  ggsave(paste0("../../out/R_analysis/plot/UMAPCluster_resolutions_",id_sample,".png"),width = 25,height = 15)
  
  # plot PCA per cell cycle -------------------------------------------------
  DimPlot(x,label = F,reduction = "pca",group.by = "Phase")
  ggsave(paste0("../../out/R_analysis/plot/PCA_",id_sample,".pdf"),width = 6,height = 4)
  ggsave(paste0("../../out/R_analysis/plot/PCA_",id_sample,".png"),width = 6,height = 4)
  
  # Elbow plot --------------------------------------------------------------
  ElbowPlot(object = x)
  ggsave(paste0("../../out/R_analysis/plot/ElbowPlot_",id_sample,".pdf"),width = 6,height = 4)
  ggsave(paste0("../../out/R_analysis/plot/ElbowPlot_",id_sample,".png"),width = 6,height = 4)
  
  # plot the high varaible genes --------------------------------------------
  plotHVG <- VariableFeaturePlot(x)
  topHVG <- head(VariableFeatures(x), 20)
  LabelPoints(plot = plotHVG, points = topHVG, repel = TRUE)
  ggsave(paste0("../../out/R_analysis/plot/HVG_",id_sample,".pdf"), width = 10, height = 6)
  ggsave(paste0("../../out/R_analysis/plot/HVG_",id_sample,".png"), width = 10, height = 6)
  
  # plot UMAPs with the clusters --------------------------------------------
  DimPlot(x,label = T)
  ggsave(paste0("../../out/R_analysis/plot/UMAP_",id_sample,".png"),width = 12,height = 9)
  ggsave(paste0("../../out/R_analysis/plot/UMAP_",id_sample,".pdf"),width = 12,height = 9)
  
  # plot UMAPs with the automatic annot -------------------------------------
  DimPlot(x,label = T,group.by = "coarse_cell_type")
  ggsave(paste0("../../out/R_analysis/plot/UMAP_autoAnnoCoarse_",id_sample,".png"),width = 12,height = 9)
  ggsave(paste0("../../out/R_analysis/plot/UMAP_autoAnnoCoarse_",id_sample,".pdf"),width = 12,height = 9)
  
  DimPlot(x,label = T,group.by = "fine_cell_type")
  ggsave(paste0("../../out/R_analysis/plot/UMAP_autoAnnoFine_",id_sample,".png"),width = 12,height = 9)
  ggsave(paste0("../../out/R_analysis/plot/UMAP_autoAnnoFine_",id_sample,".pdf"),width = 12,height = 9)
  
  # plot the number of features and the percent of mt reads -----------------
  FeaturePlot(x, features = "nFeature_RNA",cols = c("lightgrey", "red"), order = T) |
    FeaturePlot(x, features = "percent.mt",cols = c("lightgrey", "red"), order = T)
  
  ggsave(paste0("../../out/R_analysis/plot/UMAP_nFeatMito_",id_sample,".pdf"),width = 10,height = 4)
  ggsave(paste0("../../out/R_analysis/plot/UMAP_nFeatMito_",id_sample,".png"),width = 10,height = 4)
  
  # marker genes ------------------------------------------------------------
  cluster.markers <- FindAllMarkers(x,
                                    thresh.use = 0.25,
                                    test.use="wilcox",
                                    min.pct=0.25,
                                    min.diff.pct=-Inf,
                                    only.pos=TRUE)
  # save the table of markers
  write_tsv(cluster.markers,paste0("../../out/R_analysis/table/FindAllMarkers_",id_sample,".tsv"))
  
  top10 <- cluster.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  
  # SOs = ScaleData(SO, features = top10$gene)
  
  colsCLU <- createPalette(length(levels(Idents(x))), c("#ff0000", "#00ff00", "#0000ff"))
  names(colsCLU) <- levels(Idents(x))
  
  hm <- DoHeatmap(x, 
                  features = top10$gene, 
                  group.colors = colsCLU,
                  size = 4,
                  disp.min = -1.5,
                  disp.max = 1.5,
                  raster = TRUE,
                  angle = 90) +
    scale_fill_gradientn(colours = coolwarm(200)) 
  
  
  ggsave(filename = paste0("../../out/R_analysis/plot/Heatmap_FindMarkers_",id_sample,".pdf"),hm,height = 19, width = 18)
  ggsave(filename = paste0("../../out/R_analysis/plot/Heatmap_FindMarkers_",id_sample,".png"),hm,height = 19, width = 18)
})
