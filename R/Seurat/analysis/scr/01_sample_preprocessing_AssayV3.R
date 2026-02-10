# AIM ---------------------------------------------------------------------
# Apply the filter stated in the label of the file.
# This is either a testing of the filtering results, or the actual filtering choice to go forward

# libraries ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(DoubletFinder)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v3")

# define the filtering parameters -----------------------------------------
featureLow_thr <- 1000
featureHigh_thr <- 6000
mito_thr <- 15
label <- "01000_06000_15_V3"

# read in the data --------------------------------------------------------
# location of all the raw matrices
id_sample <- dir("../data/test_introns/SoupX_default_cellranger7/")

# load the LUT of the dataset
LUT <- read_csv("../data/test_introns/LUT_samples.csv")

# load the LUT of the doublet rate estimate
df_doublet <- read_csv("../data/dublets_rate_2023.csv")

# run the processing ------------------------------------------------------
# do the preprocessing over all the dataset and save the objects
list_datasc <- lapply(id_sample,function(x){
  # to track the processing of the progress of the lapply
  print(x)
  
  # read in the matrix
  data <- Read10X(data.dir = paste0("../data/test_introns/SoupX_default_cellranger7/",x))
  
  # crete the object
  datasc <- CreateSeuratObject(counts = data, project = LUT %>%
                                 filter(sample == x) %>%
                                 pull(sample),
                               # remove low expressed features
                               # in this case keep all the features before the integration
                               min.cells = 0,
                               # remove low content barcodes
                               min.features = 200)
  
  # add the metadata
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
  datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  # add also the percentage of globin. in this dataset it is not meaningful as there is no blood
  datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
  
  # label the cells based on the mt reads content
  datasc$mt_bin <- datasc@meta.data %>%
    mutate(test = case_when(percent.mt < 1~"low",
                            percent.mt < 10~"mid",
                            T ~ "high")) %>%
    pull(test)
  
  datasc$treat <- LUT %>%
    filter(sample == x) %>%
    pull(treat)
  
  datasc$test_intron <- LUT %>%
    filter(sample == x) %>%
    pull(test_intron)
  
  datasc$doxy <- LUT %>%
    filter(sample == x) %>%
    pull(doxy)
  
  datasc$exposure <- LUT %>%
    filter(sample == x) %>%
    pull(exposure)
  
  datasc$ID <- LUT %>%
    filter(sample == x) %>%
    pull(ID)
  
  # add the filtering variable based on the fixed threshold
  datasc$test <- datasc@meta.data %>%
    mutate(test = percent.mt < mito_thr & nFeature_RNA > featureLow_thr & nFeature_RNA < featureHigh_thr) %>% 
    pull(test)
  
  # add the filtering variable based on the
  stats <- cbind(log10(datasc@meta.data$nCount_RNA), log10(datasc@meta.data$nFeature_RNA),
                 datasc@meta.data$percent.mt)
  
  # add the filtering variable based on the adaptive threshold
  # library(robustbase)
  # library(scater)
  stats <- cbind(log10(datasc@meta.data$nCount_RNA),
                 log10(datasc@meta.data$nFeature_RNA),
                 datasc@meta.data$percent.mt)
  
  outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
  multi.outlier <- isOutlier(outlying, type = "higher")
  
  datasc$not_outlier <- !as.vector(multi.outlier)
  
  return(datasc)
}) %>%
  setNames(id_sample)

# confirm the class of the objects ----------------------------------------
lapply(list_datasc, function(x){
  class(x@assays$RNA)
})


# save the full metadata --------------------------------------------------
meta_total <- lapply(list_datasc, function(x){
  x@meta.data %>%
    rownames_to_column("barcode") %>%
    mutate(barcode = paste0(barcode,"|",orig.ident))
}) %>%
  bind_rows(.id = "dataset")

# count the cells that will pass the filter
meta_total %>%
  dplyr::count(ID)

# perform the filtering of the dataset based on thresholds ----------------
# perform the filtering based on the fixed threshold
list_datasc_fixed <- lapply(list_datasc,function(x){
  datasc_filter <- subset(x, subset = test == 1)
})


# pre-process the data for doublet detection ------------------------------
# This step is needed for DoubletFinder. For the doublet removal I need to run the full pre-processing steps.
list_datasc_fixed_norm <- lapply(list_datasc_fixed, function(x){
  datasc_filter <- x %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters() %>%
    RunUMAP(dims = 1:30) 
  
  return(datasc_filter)
})

# DoubletFinder -----------------------------------------------------------

# filter the dataset ======================================================
# determine the number of cell recovered
# df_doublet <- read_csv("data/doublets.csv")
# notice that I have update the file with doublets rate
# df_doublet <- read_csv("../../data/dublets_rate_2023.csv")

# use this step to remove potentially problematic samples
# list_datasc_fixed_norm_fix <- list_datasc_fixed_norm[!names(list_datasc_fixed_norm)%in%c("06_cr_61")]
list_datasc_fixed_norm_fix <- list_datasc_fixed_norm


# define the specific doublet rate per dataset =============================
# this is based on the table shared from 10x
# notice that in this case we should evaluate carefully how many doublets to assign to sample prepared with more than 30k cells. In this cases we are using the same rate applied to 30k
list_nExp <- lapply(list_datasc_fixed_norm_fix,function(x){
  recover_cell <- dim(x)[2]
  recover_cell
  # pick the rate with the lowest distance form the reference
  df_delta <- df_doublet %>% 
    mutate(abs_delta = abs(CellRecovered-recover_cell)) %>% 
    arrange(abs_delta)
  
  df_delta$MultipletRate[1]
})

# run the simulation for doublet finder ===================================
list_datasc_fixed_norm_doublet <- pmap(list(list_datasc_fixed_norm_fix,list_nExp),function(x,y){
  
  # pK Identification (no ground-truth) #####################################
  # older version of the doublet finder
  # sweep.res.list <- paramSweep_v3(x, PCs = 1:30, sct = FALSE) 
  sweep.res.list <- paramSweep(x, PCs = 1:30, sct = FALSE) 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats) 
  # plot the results to justify the pK choice 
  # bcmvn %>% 
  #   ggplot(aes(pK,BCmetric,group=1)) + 
  #   geom_point()+ 
  #   geom_line() 
  # save the max BCmetrics 
  pK <- bcmvn %>% 
    filter(BCmetric == max(BCmetric)) %>% 
    pull(pK) %>% 
    as.character() %>% 
    as.numeric() 
  
  # homotipic doublet proportion estimate ################################### 
  annotation <- x$seurat_clusters 
  homotipic.prop <- modelHomotypic(annotation) 
  
  # y is the expected nbumber of doublets based on the number of cells recovered 
  nExp.poi <- round(y*nrow(x@meta.data)) 
  nExp.poi.adj <- round(nExp.poi*(1-homotipic.prop)) 
  
  # run the DoubletFinder ################################################## 
  # older version of the doublet finder
  # pbmc.seurat.filteres <- doubletFinder_v3(x,PCs = 1:30,pN = 0.25,pK = pK,nExp = nExp.poi.adj,reuse.pANN = F,sct = F)
  pbmc.seurat.filteres <- doubletFinder(x,PCs = 1:30,pN = 0.25,pK = pK,nExp = nExp.poi.adj,reuse.pANN = F,sct = F)
  
  return(pbmc.seurat.filteres)
})

# trimm the datasets and save the filtered objects ------------------------
# fix the column name to homogenize the metadata from different objects
old_colnames <- colnames(list_datasc_fixed_norm_doublet[[1]]@meta.data)
new_colnames <- data.frame(old_colnames) %>%
  mutate(new_colnames = case_when(str_detect(old_colnames,"^pANN")~"pANN",
                                  str_detect(old_colnames,"^DF.classification")~"DF",
                                  T~old_colnames)) %>%
  pull(new_colnames)

# QC count the cells per object
lapply(list_datasc_fixed_norm_doublet, function(x){
  # check the colnames
  print(head(x@meta.data))
  # count the number of cells
  dim(x@meta.data)[1]
})

# before doublet removal ==================================================
# save the list of objects before doublet removal
list_norm_doublet <- lapply(list_datasc_fixed_norm_doublet,function(x){
  # fix the meta and save the object
  meta <- x@meta.data
  colnames(meta) <- new_colnames
  x@meta.data <- meta
  
  return(x)
}) %>% 
  setNames(names(list_datasc_fixed_norm_doublet))

saveRDS(object = list_norm_doublet,file = paste0("../out/test_introns/object/list_datasc_fix_filter_norm_doublet_SoupX_",label,".rds"))

# save the total meta with the doublet imputation
meta_total_doublet <- lapply(list_datasc_fixed_norm_doublet, function(x){
  meta <- x@meta.data
  colnames(meta) <- new_colnames
  x@meta.data <- meta
  
  x@meta.data
}) %>%
  bind_rows(.id = "dataset") %>%
  rownames_to_column("barcode")

# save the total meta
meta_total_doublet %>%
  write_tsv(paste0("../out/test_introns/table/meta_datasc_fix_filter_norm_doublet_SoupX_",label,".tsv"))

# QC count the cells per object
lapply(list_norm_doublet, function(x){
  # check the colnames
  print(head(x@meta.data))
  # count the number of cells
  dim(x@meta.data)[1]
})

# after doublet removal ===================================================
# subset the datasets after removal of the doublets before the integration
list_norm_singlets <- lapply(list_datasc_fixed_norm_doublet,function(x){
  # fix the meta
  meta <- x@meta.data
  colnames(meta) <- new_colnames
  x@meta.data <- meta
  # filter only the singlets
  datasc_filter <- subset(x, subset = DF == "Singlet")
  
  return(datasc_filter)
}) %>%
  setNames(names(list_datasc_fixed_norm_doublet))

saveRDS(object = list_norm_singlets,file = paste0("../out/test_introns/object/list_datasc_fix_filter_norm_doubletSinglet_SoupX_",label,".rds"))

# QC count the cells per object
lapply(list_norm_singlets, function(x){
  # check the colnames
  print(head(x@meta.data))
  # count the number of cells
  dim(x@meta.data)[1]
})
