# AIM ---------------------------------------------------------------------
# generate plots for the defined object

# libraries ---------------------------------------------------------------
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

# renv integration --------------------------------------------------------

# to load the packages
source(".Rprofile")

# in the config specify the following
# renv_library_path: "renv/library/linux-rocky-9.5/R-4.5/x86_64-conda-linux-gnu"

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
# read in the integrated object
data.combined <- readRDS("out/object/analysis_R45/04_sobj_filtered_harmony_afterQC.rds")

# define the label for the image names
label <- "afterQC"

# define the reference resolution
col_meta <- "RNA_snn_res.0.4"

# wrangling ---------------------------------------------------------------
# color for the full palette

color_id <- alphabet(length(unique(data.combined@meta.data[[col_meta]])))
names(color_id) <- unique(data.combined@meta.data[[col_meta]] %>% unlist())

# color for the robust palette
# color_id2 <- alphabet(length(unique(data.combined$robust_score_subclass)))
# names(color_id2) <- unique(data.combined$robust_score_subclass %>% unlist())

show_col(color_id)
# show_col(color_id2)

p01 <- clustree(data.combined@meta.data[,grep("RNA_snn_res",colnames(data.combined@meta.data))],
                prefix = "RNA_snn_res.")
ggsave(plot = p01,filename = 
         paste0("out/plot/analysis_R45/06_UMAPCluster_tree_",label,"_V5.pdf"),width = 10,height = 10)

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

p02 <- wrap_plots(list_plot)
ggsave(plot = p02,
       filename = paste0("out/plot/analysis_R45/06_UMAPCluster_resolutions_",label,"_V5.pdf"),width = 25,height = 15)

# # try to use the information from the robust annotation of Azimuth to decide on a specific resolution
# df_meta_test <- data.combined@meta.data
# 
# # save the general annotation proposed per cluster
# df_meta_test %>%
#   group_by(RNA_snn_res.0.3,custom_cell_code_general) %>%
#   summarise(.groups = "drop") %>%
#   write_tsv("../../out/table/analysis_R45/06_res03_custom_cell_code_general.tsv")

# based on this report I decided to pick resolution 0.3
Idents(data.combined) <- col_meta

# main umap
plot03 <- DimPlot(data.combined, reduction = "umap", group.by = col_meta,label = T,raster = T)
ggsave(plot = plot03,
       paste0("out/plot/analysis_R45/06_UMAPCluster_",col_meta,"_",label,"_V5.pdf"),width = 6,height = 5)

plot03a <- DimPlot(data.combined, reduction = "umap", group.by = col_meta,label = T,raster = T,split.by = "orig.ident",ncol = 3)
ggsave(plot = plot03a,
       paste0("out/plot/analysis_R45/06_UMAPClusterSplit_",col_meta,"_",label,"_V5.pdf"),width = 13,height = 10)

# plot03c <- DimPlot(data.combined, reduction = "umap", group.by = "custom_cell_code",label = T,raster = T,split.by = "orig.ident",ncol = 3)
# plot03c@data$orig.ident %>% is.na() %>% sum()
# 
# plot03c@data %>%
#   ggplot(aes(x=umap_1,y=umap_2,col=custom_cell_code)) +
#   geom_point(size = 0.1) +
#   theme_cowplot() +
#   facet_wrap(~orig.ident) +
#   theme(strip.background = element_blank()) +
#   guides(colour = guide_legend(override.aes = list(size=3)))
# ggsave("../../out/plot/analysis_R45/06_UMAPClusterSplit_customAnno_00500_05000_10_V5.pdf",width = 13,height = 10)

# main umap cell cycle
plot04 <- DimPlot(data.combined, reduction = "umap", group.by = "Phase",raster = T)
ggsave(plot=plot04,
       paste0("out/plot/analysis_R45/06_UMAPPhase_",col_meta,"_",label,"_V5.pdf"),width = 5,height = 4)

# split by sample
df_meta <- data.combined@meta.data %>%
  rownames_to_column("barcode")

df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")

data2 <- left_join(df_UMAP,df_meta,"barcode")
data2_avg <- data2 %>% group_by(.data[[col_meta]]) %>% dplyr::select(umap_1, umap_2) %>% summarize_all(mean)

# single plot
plot051 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  arrange(percent.ribo) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=umap_1,y=umap_2,col=percent.ribo),alpha=0.5,size=0.1)+
  # geom_text_repel(data = data2_avg,aes(x = umap_1, y = umap_2,label = .data[[col_meta]])) +
  scale_color_viridis_c(option = "turbo")+
  # theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave(plot = plot051,
       filename = paste0("out/plot/analysis_R45/06_UMAPRibo_",label,"_V5.pdf"),width = 6,height = 5)

# single plot
plot052 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  arrange(percent.mt) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=umap_1,y=umap_2,col=percent.mt),alpha=0.5,size=0.1)+
  # geom_text_repel(data = data2_avg,aes(x = umap_1, y = umap_2,label = .data[[col_meta]])) +
  scale_color_viridis_c(option = "turbo")+
  # theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave(plot = plot052,
       filename = paste0("out/plot/analysis_R45/06_UMAPMito_",label,"_V5.pdf"),width = 6,height = 5)

# single plot
plot053 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  arrange(nCount_RNA) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=umap_1,y=umap_2,col=nCount_RNA),alpha=0.5,size=0.1)+
  # geom_text_repel(data = data2_avg,aes(x = umap_1, y = umap_2,label = .data[[col_meta]])) +
  scale_color_viridis_c(option = "turbo")+
  # theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave(plot = plot053,
       filename = paste0("out/plot/analysis_R45/06_UMAPnCount_",label,"_V5.pdf"),width = 6,height = 5)

# split by sample
plot06 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  mutate(res_target = factor(.data[[col_meta]])) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=umap_1,y=umap_2,col=res_target),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = umap_1, y = umap_2,label = .data[[col_meta]])) +
  facet_wrap(~orig.ident,ncol=3)+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  # theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave(plot=plot06,
       filename = paste0("out/plot/analysis_R45/06_umapClusterSplit_",col_meta,"_",label,"_V5.pdf"),width = 13,height = 10)

# proportion of cell per cluster
df_summary <- df_meta %>%
  group_by(orig.ident,pathology,.data[[col_meta]]) %>%
  summarise(n = n()) %>%
  group_by(orig.ident) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)
# write_tsv(df_summary,"../out/table/06_summary_harmonySkipIntegration_AllSoupX_00500_05000_10_V5.tsv")

# define a convenient palette of colors
# show_col(hue_pal()(7))
# RColorBrewer::display.brewer.all()
# col_pal <- RColorBrewer::brewer.pal(name = "Paired",n = 3)
# col_pal <- c("#E6E6E6","#ffff09","#c1ce08","#446d05","#053c03","#4D4D4D","#06a8ce","#033b6d","#ff0ed7","#9a0404")
# show_col(col_pal)

plot07 <- df_summary %>%
  mutate(res_target = factor(.data[[col_meta]])) %>%
  # mutate(ID = factor(ID)) %>%
  # mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  ggplot(aes(x=res_target,y=prop,fill=orig.ident))+geom_col()+theme_bw()
# scale_fill_manual(values = col_pal)
# ggsave(plot=plot07,"../out/test_introns/plot/summary_harmonySkipIntegration_AllSoupX_00500_05000_10_V5.pdf",width = 7,height = 4)

plot08 <- df_summary %>%
  mutate(res_target = factor(.data[[col_meta]])) %>%
  # mutate(ID = factor(ID)) %>%
  ggplot(aes(x=res_target,y=prop,fill=orig.ident))+geom_col(position = "dodge")+theme_bw()
# scale_fill_manual(values = col_pal)
# ggsave(plot=plot08,"../out/test_introns/plot/summaryDodge_harmonySkipIntegration_AllSoupX_00500_05000_10_V5.pdf",width = 7,height = 4)

plot09 <- df_summary %>%
  mutate(res_target = factor(.data[[col_meta]])) %>%
  # mutate(ID = factor(ID)) %>%
  ggplot(aes(x=res_target,y=prop,fill=pathology))+geom_boxplot(position = "dodge")+theme_bw()+
  scale_y_continuous(trans = "sqrt")
# scale_fill_manual(values = col_pal)
# ggsave(plot=plot09,"../out/plot/06_summaryDodgeTreat_harmonySkipIntegration_AllSoupX_00500_05000_10_V5.pdf",width = 7,height = 4)

# -------------------------------------------------------------------------
# # dotplots
# # martina asked to remove FTL
# shortlist_features_list2 <- list(
#   IMMUNE = c("AIF1","TYROBP","HLA-DRA","TREM2","CX3CR1","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
#   OLIGOLINEAGE = c("CSPG4","OLIG1","OLIG2", "PDGFRA", "SOX6", "PLP1","SOX10","BCAS1","MBP","MAG"),
#   ASTRO = c("AQP4", "GFAP", "VIM","SLC1A2","S100B"),
#   # Neu = c("SYT1")
#   NEURONS = c("GAD2", "TLE4", "CUX2","SYP", "NEFL","SYT1"),
#   NPC = c("NES", "PAX6", "SOX1"),
#   CYCLING = c("TOP2A", "CDK1", "CENPF")
# )
# 
# shortlist_features_list_long <- list(
#   IMMUNE = c("IGHG1","CD38","CD8A","CD2","SKAP1","LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
#   OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10","BCAS1","MBP","MAG"),
#   ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
#   NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
#   NPC = c("NES", "PAX6", "SOX1"),
#   CYCLING = c("TOP2A", "CDK1", "CENPF"),
#   ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
#   PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1"))

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
# Idents(data.combined) <- "RNA_snn_res.0.3"
# test_short2 <- DotPlot(data.combined, features = shortlist_features_list2, dot.scale = 8,cluster.idents = T) +
#   RotatedAxis()
# ggsave(plot = test_short2,"../out/test_introns/plot/Dotplot2_harmonySkipIntegration_AllSoupX_00500_05000_10_V5.pdf",width = 13,height = 5)
# 
# test_long <- DotPlot(data.combined, features = shortlist_features_list_long, dot.scale = 8,cluster.idents = T) +
#   RotatedAxis()
# ggsave(plot=test_long,"../out/test_introns/plot/DotplotLong_harmonySkipIntegration_AllSoupX_00500_05000_10_V5.pdf",width = 25,height = 5)
# 
# df_meta %>%
#   group_by(ID) %>%
#   summarise(n=n())

# Identify conserved cell type markers ------------------------------------
# data
DefaultAssay(data.combined) <- "RNA"
Idents(data.combined) <- col_meta

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv(paste0("out/table/analysis_R45/06_FindAllMarkers_",col_meta,"_",label,"_V5.tsv"))

# pick the top 100 markers per cluster
sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  write_tsv(paste0("out/table/analysis_R45/06_FindAllMarkers_",col_meta,"_",label,"_top100_V5.tsv"))

# remove the technical genes
sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  write_tsv(paste0("out/table/analysis_R45/06_FindAllMarkers_",col_meta,"_",label,"_top100_noRIBOandMT_V5.tsv"))

# try plotting the top markers
top_specific_markers <- sobj_total_h.markers %>%
  # filter ribosomal and mt genes
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) %>%
  mutate(cluster = as.numeric(as.character(cluster))) %>%
  arrange(cluster)

# And generate e.g. a dotplot:
dittoSeq::dittoDotPlot(data.combined,
                       vars = unique(top_specific_markers$gene), 
                       group.by = col_meta) +
  scale_color_viridis_c(option = "turbo",name="relative \nexpression")
ggsave(paste0("out/plot/analysis_R45/06_Ditto_",col_meta,"_",label,"_V5.pdf"),width = 15,height = 5)
