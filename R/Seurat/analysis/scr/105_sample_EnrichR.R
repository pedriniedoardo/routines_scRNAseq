# AIM ---------------------------------------------------------------------
# run EnrichR on the table of DE from sc analysis

# librarues ---------------------------------------------------------------
library(tidyverse)
library(enrichR)
library(scales)

# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
#filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Azimuth"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_Pathways_2024")

# split by direction ------------------------------------------------------

res_file <- "../out/table/100_response_STIM_vs_CTRL_seurat_sc_V5.tsv"

list_genes_UP <- read_tsv(res_file) %>%
  filter(p_val_adj<0.05 & avg_log2FC > 1&!is.na(gene)) %>%
  split(f = .$annotation) %>%
  map(function(x){
    x %>%
      pull(gene)
  })

list_genes_DOWN <- read_tsv(res_file) %>%
  filter(p_val_adj<0.05 & avg_log2FC < (-1)&!is.na(gene)) %>%
  split(f = .$annotation) %>%
  map(function(x){
    x %>%
      pull(gene)
  })

# define the background
# background <- df_modules$feature

# x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
list_enrichr_UP <- lapply(list_genes_UP,function(x){
  genes <- x
  # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>% 
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_enrichr_UP %>%
  write_tsv("../out/table/105_enrichR_DE_UP.tsv")

list_enrichr_DOWN <- lapply(list_genes_DOWN,function(x){
  genes <- x
  # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>% 
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_enrichr_DOWN %>%
  write_tsv("../out/table/105_enrichR_DE_DOWN.tsv")

# build the plots
plot_list_UP <- list_enrichr_UP %>%
  split(f = .$comparison)

plot_list_DOWN <- list_enrichr_DOWN %>%
  split(f = .$comparison)

# library(scales)
list_plot_UP <- pmap(list(plot_list_UP,names(plot_list_UP)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")
  
  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_UP)
ggsave("../out/plot/105_enrichR_DE_UP.pdf",width = 36,height = 36,limitsize = FALSE)

list_plot_DOWN <- pmap(list(plot_list_DOWN,names(plot_list_DOWN)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")
  
  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_DOWN)
ggsave("../out/plot/105_enrichR_DE_DOWN.pdf",width = 36,height = 36,limitsize = FALSE)

# # analysis regardless of the direction ------------------------------------
# # list of genes to consider for the enrichment analysis
# test_out <- read_tsv("../out/table/100_response_STIM_vs_CTRL_seurat_sc_V5.tsv") %>%
#   # filter(cluster %in% "CD14 Mono") %>%
#   filter(p_val_adj < 0.05) %>%
#   filter(abs(avg_log2FC) > 1) %>%
#   split(f = .$annotation)
# glimpse(test_out)
# 
# # read in the metadata for the genes
# # meta <- read_csv("data/scTrem/GSE130626_gene_info.csv")
# 
# #
# list_res_tot <- lapply(test_out, function(x){
#   x %>%
#     pull(gene)
# })
# 
# # query -------------------------------------------------------------------
# list <- lapply(list_res_tot,function(x){
#   genes <- x
#   out_enrich <- enrichr(genes, dbs_db)
#   #
#   out_enrich %>%
#     bind_rows(.id = "annotation")
# })
# 
# df_enrichr_annotation_enriched_tot <- list %>%
#   bind_rows(.id = "comparison")
# 
# df_enrichr_annotation_enriched_tot %>%
#   write_tsv("../out/table/03_enrichR_deseq2.tsv")
# 
# # check entries for senescence
# df_enrichr_annotation_enriched_tot %>%
#   filter(str_detect(Term,pattern = "enescence"))
# 
# # library(scales)
# list_plot <- lapply(unique(df_enrichr_annotation_enriched_tot$comparison),function(x){
#   df_enrichr_annotation_enriched_tot %>%
#     filter(comparison == x) %>%
#     group_by(annotation) %>%
#     arrange(P.value) %>%
#     dplyr::slice(1:20) %>%
#     mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
#     mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
#     ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
#     theme(strip.background = element_blank()) +
#     scale_color_gradientn(colors = c("red","blue"),
#                           values = rescale(c(0,1)),
#                           limits = c(0,0.2))
#   
# }) %>%
#   setNames(unique(df_enrichr_annotation_enriched_tot$comparison))
# 
# # scale_color_gradient(low = "red",high = "blue")
# pmap(list(list_plot,names(list_plot)), function(x,y){
#   ggsave(plot = x,paste0("../out/plot/03_enrichR_deseq2_",y,".pdf"),width = 7,height = 10)
# })
