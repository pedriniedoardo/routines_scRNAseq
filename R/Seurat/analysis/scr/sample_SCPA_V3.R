# libraries ---------------------------------------------------------------
library(tidyverse)
library(SCPA)
library(msigdbr)
library(Seurat)
library(ComplexHeatmap)
library(circlize)

# I/O ---------------------------------------------------------------------
input <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(input,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(input,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(input,label = T,raster = T,group.by = "treat_full")

# wrangling ---------------------------------------------------------------
# add the treatmend varibale to the metadata
meta <- input@meta.data %>%
  rownames_to_column("barcodes")

# meta_full <- left_join(meta,LUT,by=c("orig.ident"="sample"))
# # add tot he
# input$treat <- meta_full$disease_fix

# explore the metadata
table(input@meta.data$treat_full)
table(input@meta.data$expertAnno.l1)
table(input@meta.data$harmonized_donor2)
table(input@meta.data$expertAnno.l1,input@meta.data$harmonized_donor2)

# create the splitting variable for the dataset if not already present

# split the seurat object into the different cell types
Idents(input) <- "expertAnno.l1"
list_input <- SplitObject(input,split.by = "expertAnno.l1")

# calculation -------------------------------------------------------------
# define the annotation
# pathways <- c("hallmark", "kegg", "reactome")
pathways_id <- c("kegg")
pathways <- msigdbr("Homo sapiens") %>%
  filter(grepl(paste(pathways, collapse = "|"), gs_name, ignore.case = T)) %>%
  format_pathways()

# Id for looping over the cell type
id_cellID <- levels(Idents(input))
# Id for looping over the treatment condition in this case use BASELINE as reference
id_treat <- c("CSF.ctrl.24h","CSF.MS.24h","CSF.MS.48h","cytokine","Fe","myelin","TBHP")

# for testing
# cellID <- "MG"
# treat <- "Fe"

# run the comparison all vs BASELINE pairwise
SCPA_BASELINE <- lapply(id_cellID,function(cellID){
  # keep track of the progress
  print(cellID)
  df_cellID <- lapply(id_treat,function(treat){
    # keep track of the progress
    print(paste(treat,cellID))
    
    # define the input dataset
    data <- list_input[[cellID]]
    
    # define the comparison
    control <- seurat_extract(data,meta1 = "treat_full", value_meta1 = "BASELINE")
    disease <- seurat_extract(data,meta1 = "treat_full", value_meta1 = treat)
    
    # run the comparison
    scpa_out <- compare_pathways(list(control, disease), pathways,parallel = T,cores = 8)
    scpa_out %>%
      mutate(comparison = paste0(treat,"_vs_BASELINE")) %>%
      mutate(cellID = cellID)
    # save the output
    # write_tsv(x = scpa_out,file = paste0("out/table/SCPA_out/SCPA_Ms_vs_Ctrl_",x,"_Jakel.tsv"))
    
  }) %>%
    bind_rows()
  # save the table version
  write_tsv(x = df_cellID,file = paste0("../../out/table/SCPA_out/SCPA_vs_BASELINE_",cellID,".tsv"))
  return(df_cellID)
}) %>%
  bind_rows()

write_tsv(x = SCPA_BASELINE,file = "../../out/table/SCPA_out/SCPA_vs_BASELINE_all.tsv")

# do the same but use as reference the CSF ctrl
# Id for looping over the treatment condition in this case use BASELINE as reference
id_treat2 <- c("BASELINE","CSF.MS.24h","CSF.MS.48h","cytokine","Fe","myelin","TBHP")

# run the comparison all vs CSFctrl pairwise
SCPA_CSFctrl <- lapply(id_cellID,function(cellID){
  # keep track of the progress
  print(cellID)
  df_cellID <- lapply(id_treat2,function(treat){
    # keep track of the progress
    print(paste(treat,cellID))
    
    # define the input dataset
    data <- list_input[[cellID]]
    
    # define the comparison
    control <- seurat_extract(data,meta1 = "treat_full", value_meta1 = "CSF.ctrl.24h")
    disease <- seurat_extract(data,meta1 = "treat_full", value_meta1 = treat)
    
    # run the comparison
    scpa_out <- compare_pathways(list(control, disease), pathways,parallel = T,cores = 8)
    scpa_out %>%
      mutate(comparison = paste0(treat,"_vs_CSFctrl")) %>%
      mutate(cellID = cellID)
    # save the output
    # write_tsv(x = scpa_out,file = paste0("out/table/SCPA_out/SCPA_Ms_vs_Ctrl_",x,"_Jakel.tsv"))
    
  }) %>%
    bind_rows()
  # save the table version
  write_tsv(x = df_cellID,file = paste0("../../out/table/SCPA_out/SCPA_vs_CSFctrl_",cellID,".tsv"))
  return(df_cellID)
}) %>%
  bind_rows()

write_tsv(x = SCPA_CSFctrl,file = "../../out/table/SCPA_out/SCPA_vs_CSFctrl_all.tsv")

# Plotting a global summary of the data -----------------------------------
# pull the file with all the cell types
# If you're doing a 2 sample comparison like your example above, there will be a fold change column in the output from compare_pathways, which will tell you the directionality. The fold change is calculated from population1 - population2, so a negative value will be higher in population2
SCPA_BASELINE <- read_tsv(file = "../../out/table/SCPA_out/SCPA_vsBESELINE_all.tsv")
SCPA_CSFctrl <- read_tsv(file = "../../out/table/SCPA_out/SCPA_vs_CSFctrl_all.tsv")

# plotting ----------------------------------------------------------------

# run the plotting with BASELINE as reference -----------------------------
# default scatter plot
SCPA_BASELINE <- SCPA_BASELINE %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

# pull the lebels. use the top 5 per comparision
SCPA_BASELINE_label <- SCPA_BASELINE %>%
  filter(color == "mediumseagreen") %>%
  group_by(comparison,cellID) %>%
  arrange(comparison,cellID,desc(qval)) %>%
  dplyr::slice(1:5) %>%
  # shorten the label
  mutate(pathway2 = str_remove(Pathway,pattern = "KEGG_"))

# aa_path <- rest_act %>% 
#   filter(grepl(pattern = "reactome_arachi", ignore.case = T, x = Pathway))
SCPA_BASELINE %>%
  ggplot(aes(-FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = SCPA_BASELINE$color, stroke = 0.3) +
  # geom_point(data = aa_path, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
  # xlim(-20, 80) +
  ylim(0, 15) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)+
  facet_grid(comparison~cellID)+
  theme(strip.background = element_blank()) +
  ggrepel::geom_text_repel(data = SCPA_BASELINE_label,aes(x=-FC,y=qval,label = pathway2),force = 200,min.segment.length = 0,segment.alpha=0.5,nudge_x = 2,nudge_y = 2)
ggsave("../../out/image/SCPA_vsBESELINE_all_scatter.pdf",width = 40,height = 40)

# default heatmap
df_hm_baseline <- SCPA_BASELINE %>%
  mutate(group = paste0(comparison,"|",cellID)) %>%
  dplyr::select(Pathway,qval,group) %>%
  pivot_wider(names_from = group,values_from = qval) %>% 
  mutate(Pathway = str_remove_all(Pathway,pattern = "KEGG_")) %>%
  column_to_rownames("Pathway")

# annotation of the column
meta_baseline <- data.frame(group = colnames(df_hm_baseline)) %>%
  separate(group,into = c("treat","cellID"),sep = "\\|",remove = F)

# hue_pal()(7)
# brewer_pal(palette = "Set1")(7)
col_an_baseline <- HeatmapAnnotation(cellID = meta_baseline$cellID,
                                     comparison = meta_baseline$treat,
                                     col = list(cellID = c("ASTRO" = "#F8766D",
                                                           "GLIA_IMM" = "#C49A00",
                                                           "MG" = "#53B400",
                                                           "NEU" = "#00C094",
                                                           "OLIGO" = "#00B6EB",
                                                           "OPC" = "#A58AFF",
                                                           "PROG" = "#FB61D7"),
                                                comparison = c("CSF.ctrl.24h_vs_BASELINE" = "#FFFF33",
                                                               "CSF.MS.24h_vs_BASELINE" = "#FF7F00",
                                                               "CSF.MS.48h_vs_BASELINE" = "#E41A1C",
                                                               "cytokine_vs_BASELINE" = "#984EA3",
                                                               "Fe_vs_BASELINE" = "#A65628",
                                                               "myelin_vs_BASELINE" = "#377EB8",
                                                               "TBHP_vs_BASELINE" = "#4DAF4A"))
)



hm_baseline <- df_hm_baseline %>%
  Heatmap(name = "Qval",
          show_row_names = F, 
          top_annotation = col_an_baseline,
          border = T,
          show_row_dend = F,
          show_column_dend = F,
          show_column_names = F)
pdf("../../out/image/SCPA_vsBESELINE_all_heatmap.pdf",width = 6,height = 6) 
draw(hm_baseline,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# plot with labels on the row
hm_baseline2 <- df_hm_baseline %>%
  Heatmap(name = "Qval",
          show_row_names = T, 
          row_names_gp = gpar(fontsize = 3),
          top_annotation = col_an_baseline,
          border = T,
          show_row_dend = F,
          show_column_dend = F,
          show_column_names = F)
pdf("../../out/image/SCPA_vsBESELINE_all_heatmap2.pdf",width = 6,height = 10) 
draw(hm_baseline2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# Finding something a bit more specific
# You obviously now want to use this data to find something a bit more biologically relevant. One way to do this would be to see if any pathways were tissue specific. To do this, we calculated the most variable pathways across comparisons, and then plotted the variance against the pathway (pathway names added later). Here we can see signatures for antimicrobial peptide production and prostaglandin synthesis at highly variable across conditions.
plot_variance_baseline <- apply(df_hm_baseline, 1, var) %>%
  data.frame(variance = .) %>% 
  arrange(desc(variance)) %>% 
  rownames_to_column("pathway") %>%
  mutate(rank = nrow(.):1)

plot_variance_baseline %>%
  ggplot(aes(rank, variance)) +
  geom_point(shape = 21, cex = 3, fill = "royalblue2", color = 'black', stroke = 0.2) +
  scale_x_discrete(expand = c(0.04, 0.04)) +
  labs(x = "Pathway", y = "Variance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))+
  ggrepel::geom_text_repel(data = plot_variance_baseline %>% dplyr::slice(1:10),
                           aes(x=rank,y=variance,label = pathway),
                           force = 200,min.segment.length = 0,segment.alpha=0.5,nudge_x = 2,nudge_y = 2)
ggsave("../../out/image/SCPA_vsBESELINE_all_variance.pdf",width = 10,height = 10)

# run the same using as baseline CSF ctrl ---------------------------------
SCPA_CSFctrl <- SCPA_CSFctrl %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

# pull the lebels. use the top 5 per comparision
SCPA_CSFctrl_label <- SCPA_CSFctrl %>%
  filter(color == "mediumseagreen") %>%
  group_by(comparison,cellID) %>%
  arrange(comparison,cellID,desc(qval)) %>%
  dplyr::slice(1:5) %>%
  # shorten the label
  mutate(pathway2 = str_remove(Pathway,pattern = "KEGG_"))

# aa_path <- rest_act %>% 
#   filter(grepl(pattern = "reactome_arachi", ignore.case = T, x = Pathway))
SCPA_CSFctrl %>%
  ggplot(aes(-FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = SCPA_CSFctrl$color, stroke = 0.3) +
  # geom_point(data = aa_path, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
  # xlim(-20, 80) +
  ylim(0, 15) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)+
  facet_grid(comparison~cellID)+
  theme(strip.background = element_blank()) +
  ggrepel::geom_text_repel(data = SCPA_CSFctrl_label,aes(x=-FC,y=qval,label = pathway2),force = 200,min.segment.length = 0,segment.alpha=0.5,nudge_x = 2,nudge_y = 2)
ggsave("../../out/image/SCPA_vsCSFctrl_all_scatter.pdf",width = 40,height = 40)

# default heatmap
df_hm_csf <- SCPA_CSFctrl %>%
  mutate(group = paste0(comparison,"|",cellID)) %>%
  dplyr::select(Pathway,qval,group) %>%
  pivot_wider(names_from = group,values_from = qval) %>% 
  mutate(Pathway = str_remove_all(Pathway,pattern = "KEGG_")) %>%
  column_to_rownames("Pathway")

# annotation of the column
meta_csf <- data.frame(group = colnames(df_hm_csf)) %>%
  separate(group,into = c("treat","cellID"),sep = "\\|",remove = F)

# hue_pal()(7)
# brewer_pal(palette = "Set1")(7)
col_an_csf <- HeatmapAnnotation(cellID = meta_csf$cellID,
                                comparison = meta_csf$treat,
                                col = list(cellID = c("ASTRO" = "#F8766D",
                                                      "GLIA_IMM" = "#C49A00",
                                                      "MG" = "#53B400",
                                                      "NEU" = "#00C094",
                                                      "OLIGO" = "#00B6EB",
                                                      "OPC" = "#A58AFF",
                                                      "PROG" = "#FB61D7"),
                                           comparison = c("BASELINE_vs_CSFctrl" = "#FFFF33",
                                                          "CSF.MS.24h_vs_CSFctrl" = "#FF7F00",
                                                          "CSF.MS.48h_vs_CSFctrl" = "#E41A1C",
                                                          "cytokine_vs_CSFctrl" = "#984EA3",
                                                          "Fe_vs_CSFctrl" = "#A65628",
                                                          "myelin_vs_CSFctrl" = "#377EB8",
                                                          "TBHP_vs_CSFctrl" = "#4DAF4A"))
)



hm_csf <- df_hm_csf %>%
  Heatmap(name = "Qval",
          show_row_names = F, 
          top_annotation = col_an_csf,
          border = T,
          show_row_dend = F,
          show_column_dend = F,
          show_column_names = F)
pdf("../../out/image/SCPA_vsCSFctrl_all_heatmap.pdf",width = 6,height = 6) 
draw(hm_csf,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# plot with labels on the row
hm_csf2 <- df_hm_csf %>%
  Heatmap(name = "Qval",
          show_row_names = T, 
          row_names_gp = gpar(fontsize = 3),
          top_annotation = col_an_csf,
          border = T,
          show_row_dend = F,
          show_column_dend = F,
          show_column_names = F)
pdf("../../out/image/SCPA_vsCSFctrl_all_heatmap2.pdf",width = 6,height = 10) 
draw(hm_csf2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# Finding something a bit more specific
# You obviously now want to use this data to find something a bit more biologically relevant. One way to do this would be to see if any pathways were tissue specific. To do this, we calculated the most variable pathways across comparisons, and then plotted the variance against the pathway (pathway names added later). Here we can see signatures for antimicrobial peptide production and prostaglandin synthesis at highly variable across conditions.
plot_variance_csf <- apply(df_hm_csf, 1, var) %>%
  data.frame(variance = .) %>% 
  arrange(desc(variance)) %>% 
  rownames_to_column("pathway") %>%
  mutate(rank = nrow(.):1)

plot_variance_csf %>%
  ggplot(aes(rank, variance)) +
  geom_point(shape = 21, cex = 3, fill = "royalblue2", color = 'black', stroke = 0.2) +
  scale_x_discrete(expand = c(0.04, 0.04)) +
  labs(x = "Pathway", y = "Variance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))+
  ggrepel::geom_text_repel(data = plot_variance_csf %>% dplyr::slice(1:10),
                           aes(x=rank,y=variance,label = pathway),
                           force = 200,min.segment.length = 0,segment.alpha=0.5,nudge_x = 2,nudge_y = 2)
ggsave("../../out/image/SCPA_vsCSFctrl_all_variance.pdf",width = 10,height = 10)


