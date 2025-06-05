# AIM ---------------------------------------------------------------------
# the aim is to plot the data after integration.
# this integration was run by skipping the seurat integration (by merging the matrices) and running Harmony

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
library(patchwork)
library(Nebulosa)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
data.combined <- readRDS("../out/test_introns/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.rds")
Idents(data.combined) <- "RNA_snn_res.0.5"

# plots -------------------------------------------------------------------
# library(clustree)
clustree::clustree(data.combined@meta.data[,grep("RNA_snn_res",colnames(data.combined@meta.data))],
                   prefix = "RNA_snn_res.")
ggsave("../out/test_introns/plot/UMAPCluster_tree_V5.pdf",width = 10,height = 10)

# plot the UMAP with all the resolutions runs
id_resolution <- str_subset(colnames(data.combined@meta.data),pattern = "RNA_snn_res") %>%
  sort()

list_plot <- lapply(id_resolution,function(x){
  plot <- DimPlot(data.combined,
                  reduction = "umap",
                  group.by = x,
                  label = T,
                  raster = T)
  return(plot)
})

wrap_plots(list_plot)
ggsave("../out/test_introns/plot/UMAPCluster_resolutions_V5.pdf",width = 25,height = 15)

# main umap
plot03 <- DimPlot(data.combined, reduction = "umap", group.by = "RNA_snn_res.0.5",label = T,raster = T)
ggsave(plot = plot03,"../out/test_introns/plot/UMAPCluster_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.pdf",width = 6,height = 5)

plot03a <- DimPlot(data.combined, reduction = "umap", group.by = "RNA_snn_res.0.5",label = T,raster = T,split.by = "ID",ncol = 3)
ggsave(plot = plot03a,"../out/test_introns/plot/UMAPClusterSplit_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.pdf",width = 11,height = 5)

# main umap cell cycle
plot04 <- DimPlot(data.combined, reduction = "umap", group.by = "Phase",raster = T)
ggsave(plot=plot04,"../out/test_introns/plot/UMAPPhase_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.pdf",width = 5,height = 4)

# split by sample
df_meta <- data.combined@meta.data %>%
  rownames_to_column("barcode")

df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")

data2 <- left_join(df_UMAP,df_meta,"barcode")
data2_avg <- data2 %>% group_by(RNA_snn_res.0.5) %>% dplyr::select(umap_1, umap_2) %>% summarize_all(mean)

# single plot
# plot05 <- data2 %>%
#   mutate(orig.ident = factor(orig.ident)) %>%
#   # arrange(signature_score1) %>%
#   # mutate(gene = "Ptx3") %>%
#   ggplot() + 
#   geom_point(aes(x=umap_1,y=umap_2,col=RNA_snn_res.0.5),alpha=0.5,size=0.1)+
#   geom_text_repel(data = data2_avg,aes(x = umap_1, y = umap_2,label = RNA_snn_res.0.5)) +
#   guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
#   theme(strip.background = element_blank(), 
#         panel.border = element_rect(colour = "black", fill = NA))+
#   theme_bw()
# ggsave(plot=plot05,"../../out/plot/UMAPClusterggplot_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 6,height = 4)
# ggsave(plot=plot05,"../../out/plot/UMAPClusterggplot_harmonySkipIntegration_AllSoupX_01000_06000_15.png",width = 6,height = 4)

# save the same with theme cowplot
# plot05_alt <- data2 %>%
#   mutate(orig.ident = factor(orig.ident)) %>%
#   # arrange(signature_score1) %>%
#   # mutate(gene = "Ptx3") %>%
#   ggplot() + 
#   geom_point(aes(x=umap_1,y=umap_2,col=RNA_snn_res.0.5),alpha=0.5,size=0.1)+
#   geom_text_repel(data = data2_avg,aes(x = umap_1, y = umap_2,label = RNA_snn_res.0.5)) +
#   guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
#   theme(strip.background = element_blank(), 
#         panel.border = element_rect(colour = "black", fill = NA))+
#   theme_cowplot()
# ggsave(plot=plot05_alt,"../../out/plot/UMAPClusterggplot2_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 6,height = 4)
# ggsave(plot=plot05_alt,"../../out/plot/UMAPClusterggplot2_harmonySkipIntegration_AllSoupX_01000_06000_15.png",width = 6,height = 4)

# single plot
plot051 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  arrange(percent.ribo) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=umap_1,y=umap_2,col=percent.ribo),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = umap_1, y = umap_2,label = RNA_snn_res.0.5)) +
  scale_color_viridis_c(option = "turbo")+
  # theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave(plot = plot051,"../out/test_introns/plot/UMAPRibo_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.pdf",width = 5,height = 4)
ggsave(plot = plot051,"../out/test_introns/plot/UMAPRibo_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.png",width = 5,height = 4,bg = "white")

# single plot
plot052 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  arrange(percent.mt) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=umap_1,y=umap_2,col=percent.mt),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = umap_1, y = umap_2,label = RNA_snn_res.0.5)) +
  scale_color_viridis_c(option = "turbo")+
  # theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave(plot = plot052,"../out/test_introns/plot/UMAPMito_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.pdf",width = 5,height = 4)
ggsave(plot = plot052,"../out/test_introns/plot/UMAPMito_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.png",width = 5,height = 4,bg="white")

# single plot
plot053 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  arrange(nCount_RNA) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=umap_1,y=umap_2,col=nCount_RNA),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = umap_1, y = umap_2,label = RNA_snn_res.0.5)) +
  scale_color_viridis_c(option = "turbo")+
  # theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave(plot = plot053,"../out/test_introns/plot/UMAPnCount_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.pdf",width = 5,height = 4)
ggsave(plot = plot053,"../out/test_introns/plot/UMAPnCount_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.png",width = 5,height = 4,bg="white")

# plot densities using the nubulosa package
# deafult plotting the metadata
# plot_density(data.combined,reduction = "umap",features = "percent.ribo")
# default ploting genes
# plot_density(data.combined,reduction = "umap",features = c("GFAP","RORB","CD38"))
# plot changing the color palette
plot_density(data.combined,reduction = "umap",features = "percent.ribo")+ggplot2::scale_color_viridis_c(option = "turbo")

# split by sample
plot06 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=umap_1,y=umap_2,col=RNA_snn_res.0.5),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = umap_1, y = umap_2,label = RNA_snn_res.0.5)) +
  facet_wrap(~ID,ncol=3)+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  # theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave(plot=plot06,"../out/test_introns/plot/umapClusterSplit_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.pdf",width = 11,height = 5)
ggsave(plot=plot06,"../out/test_introns/plot/umapClusterSplit_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.png",width = 11,height = 5,bg="white")

# proportion of cell per cluster
df_summary <- df_meta %>%
  group_by(ID,RNA_snn_res.0.5) %>%
  summarise(n = n()) %>%
  group_by(ID) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)
write_tsv(df_summary,"../out/test_introns/table/summary_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.tsv")

# define a convenient palette of colors
show_col(hue_pal()(7))
RColorBrewer::display.brewer.all()
col_pal <- RColorBrewer::brewer.pal(name = "Paired",n = 3)
# col_pal <- c("#E6E6E6","#ffff09","#c1ce08","#446d05","#053c03","#4D4D4D","#06a8ce","#033b6d","#ff0ed7","#9a0404")
show_col(col_pal)

plot07 <- df_summary %>%
  mutate(ID = factor(ID)) %>%
  # mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  ggplot(aes(x=RNA_snn_res.0.5,y=prop,fill=ID))+geom_col()+theme_bw()+
  scale_fill_manual(values = col_pal)
ggsave(plot=plot07,"../out/test_introns/plot/summary_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.pdf",width = 7,height = 4)

plot08 <- df_summary %>%
  mutate(ID = factor(ID)) %>%
  ggplot(aes(x=RNA_snn_res.0.5,y=prop,fill=ID))+geom_col(position = "dodge")+theme_bw()+scale_fill_manual(values = col_pal)
ggsave(plot=plot08,"../out/test_introns/plot/summaryDodge_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.pdf",width = 7,height = 4)

# dotplots
# martina asked to remove FTL
shortlist_features_list2 <- list(
  IMMUNE = c("AIF1","TYROBP","HLA-DRA","TREM2","CX3CR1","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("CSPG4","OLIG1","OLIG2", "PDGFRA", "SOX6", "PLP1","SOX10","BCAS1","MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "VIM","SLC1A2","S100B"),
  # Neu = c("SYT1")
  NEURONS = c("GAD2", "TLE4", "CUX2","SYP", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF")
)

shortlist_features_list_long <- list(
  IMMUNE = c("IGHG1","CD38","CD8A","CD2","SKAP1","LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10","BCAS1","MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
  NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF"),
  ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
  PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1"))

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
Idents(data.combined) <- "RNA_snn_res.0.5"
test_short2 <- DotPlot(data.combined, features = shortlist_features_list2, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
ggsave(plot = test_short2,"../out/test_introns/plot/Dotplot2_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.pdf",width = 13,height = 5)

test_long <- DotPlot(data.combined, features = shortlist_features_list_long, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
ggsave(plot=test_long,"../out/test_introns/plot/DotplotLong_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.pdf",width = 25,height = 5)

df_meta %>%
  group_by(ID) %>%
  summarise(n=n())

# Identify conserved cell type markers ------------------------------------
# data
DefaultAssay(data.combined) <- "RNA"

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("../out/test_introns/table/FindAllMarkers_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.tsv")

# pick the top 100 markers per cluster
sobj_total_h.markers %>%
  group_by(cluster) %>%
  slice(1:100) %>%
  write_tsv("../out/test_introns/table/FindAllMarkers_harmonySkipIntegration_AllSoupX_01000_06000_15_top100_V5.tsv")

# remove the technical genes
sobj_total_h.markers %>%
  group_by(cluster) %>%
  slice(1:100) %>%
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  write_tsv("../out/test_introns/table/FindAllMarkers_harmonySkipIntegration_AllSoupX_01000_06000_15_top100_noRIBOandMT_V5.tsv")

# try plotting the top markers
top_specific_markers <- sobj_total_h.markers %>%
  # filter ribosomal and mt genes
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)

# And generate e.g. a dotplot:
dittoSeq::dittoDotPlot(data.combined,
                       vars = unique(top_specific_markers$gene), 
                       group.by = "RNA_snn_res.0.5")+scale_color_viridis_c(option = "turbo",name="relative \nexpression")
ggsave("../out/test_introns/plot/Ditto_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.pdf",width = 10,height = 5)
