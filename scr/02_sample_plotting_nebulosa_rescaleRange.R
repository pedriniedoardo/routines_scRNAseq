# lirbaries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ks)

# function definition -----------------------------------------------------
.extract_feature_data <- function(exp_data, features) {
  # Extract data for input features
  i <- colnames(exp_data) %in% features
  
  # Test existence of feature in gene expression data
  j <- !features %in% colnames(exp_data)
  if (any(j)) {
    stop(
      "'", paste(features[j], collapse = ", "),
      "' feature(s) not present in meta.data or expression data"
    )
  }
  vars <- exp_data[, i, drop = FALSE]
  vars <- vars[, features, drop = FALSE]
  vars
}

calculate_density <- function(w, x, method, adjust = 1, map = TRUE) {
  if (method == "ks") {
    dens <- kde(x[, c(1, 2)],
                w = w / sum(w) * length(w)
    )
  } else if (method == "wkde") {
    dens <- wkde2d(
      x = x[, 1],
      y = x[, 2],
      w = w / sum(w) * length(w),
      adjust = adjust
    )
  }
  
  if (map) {
    get_dens(x, dens, method)
  } else {
    dens
  }
}

get_dens <- function(data, dens, method) {
  if (method == "ks") {
    ix <- findInterval(data[, 1], dens$eval.points[[1]])
    iy <- findInterval(data[, 2], dens$eval.points[[2]])
    ii <- cbind(ix, iy)
    z <- dens$estimate[ii]
  } else if (method == "wkde") {
    ix <- findInterval(data[, 1], dens$x)
    iy <- findInterval(data[, 2], dens$y)
    ii <- cbind(ix, iy)
    z <- dens$z[ii]
  }
  z
}

# var definition ----------------------------------------------------------
# joint = FALSE
# dims = c(1, 2)
# method = c("ks", "wkde")
# adjust = 1
size = 1
shape = 16
# combine = TRUE
# pal = "viridis"
# raster = F

# load a sample data ------------------------------------------------------
# test <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
test <- readRDS("../out/test_introns/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15_V5.rds")

# load the siganture file
# list_sig <- readRDS("../../data/signatures/senescence_pathways.rds")

# add the module score
# x <- "senmayo"

# signature.genes.df <- list_sig[[x]]

# pull the genes
# signature.genes <- signature.genes.df %>%
#   pull(Genes) %>%
#   unique()
signature.genes <- c("AIF1","TYROBP","HLA-DRA","TREM2","CX3CR1","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC")

# score the module
test <- AddModuleScore(test,
                       features = list(signature.genes),
                       name="signature_score")

# wrangling ---------------------------------------------------------------
# split the dataset
full_data <- as.data.frame(SeuratObject::Embeddings(test[["umap"]])) %>%
  bind_cols(test@meta.data) %>%
  mutate(signature_score2 = signature_score1+abs(min(signature_score1)))
# bind_cols(SeuratObject::FetchData(test, vars = "CD4", slot = "data")) %>%

# split the dataset per batch
list_cell_embeddings <- split(full_data,f = full_data$treat) %>%
  lapply(function(x){
    x %>%
      # dplyr::select("UMAP_1","UMAP_2")
      dplyr::select("umap_1","umap_2")
  })

list_vars <- split(full_data,f = full_data$treat) %>%
  lapply(function(x){
    x %>%
      # dplyr::select("CD4")
      dplyr::select("signature_score2")
  })

# vars <- .extract_feature_data(exp_data, "CD4")
dim_names <- colnames(list_cell_embeddings$BASELINE)

# library(ks)
# measure the density splitwise
df_density <- pmap(list(list_vars,list_cell_embeddings),function(v,cell_emb){
  # vars[, 1]
  z <- calculate_density(v[, 1], cell_emb, method = "ks", adjust)
  
  # define the dataset
  df <- data.frame(cell_emb,v,density = z)
  return(df)
}) %>%
  bind_rows(.id = "treat")

head(df_density)

df_density %>%
  ggplot(aes(x=signature_score2))+geom_histogram()

df_density %>%
  ggplot(aes(x=density))+geom_histogram()

# df_density %>%
#   ggplot() +
#   aes_string(dim_names[1], dim_names[2], color = "density") +
#   geom_point(shape = shape, size = size) +
#   xlab(gsub("_", " ", dim_names[1])) +
#   ylab(gsub("_", " ", dim_names[2])) +
#   ggtitle("CD4") +
#   labs(color = "density") +
#   theme(
#     text = element_text(size = 14),
#     panel.background = element_blank(),
#     axis.text.x = element_text(color = "black"),
#     axis.text.y = element_text(color = "black"),
#     axis.line = element_line(size = 0.25),
#     strip.background = element_rect(color = "black", fill = "#ffe5cc")
#   )


# try to add a splitting variable
df_density %>%
  ggplot() +
  aes_string(dim_names[1], dim_names[2], color = "density") +
  geom_point(shape = shape, size = size) +
  xlab(gsub("_", " ", dim_names[1])) +
  ylab(gsub("_", " ", dim_names[2])) +
  ggtitle("GOI") +
  labs(color = "density") +
  theme_void() +
  theme(strip.background = element_blank()) + facet_wrap(~treat)+scale_color_viridis_c(option = "turbo")+
  scale_color_gradientn("sig score",colours = viridis::turbo(10),limits = c(0,0.06),oob = scales::squish)
# ggsave(paste0("../../out/image/revision/121_UMAP_score_",x,"_tailored2_nebulosa.pdf"),width = 11,height = 10)
# ggsave(paste0("../../out/image/revision/121_UMAP_score_",x,"_tailored2_nebulosa.png"),width = 11,height = 10,bg = "white")

# df_density %>%
#   arrange(signature_score1) %>%
#   mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Mye","CX_Demye","WM_Ctrl","WM_NAWM","WM_CA","WM_CI","WM_Core"))) %>% 
#   # mutate(gene = "Ptx3") %>%
#   ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
#   xlab(gsub("_", " ", dim_names[1])) +
#   ylab(gsub("_", " ", dim_names[2])) +
#   ggtitle("senmayo") +
#   # labs(color = "density") +
#   theme_void() +
#   theme(strip.background = element_blank()) + facet_wrap(~pathology_class)+scale_color_viridis_c(option = "turbo")+
#   scale_color_gradientn("sig score",colours = viridis::turbo(10),limits = c(-0.1,0.3),oob = scales::squish)