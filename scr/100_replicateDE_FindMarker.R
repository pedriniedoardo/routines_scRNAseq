# replicate the pvalue ----------------------------------------------------
# https://github.com/satijalab/seurat/blob/9755c164d99828dbc5dd9c8364389766cd4ff7fd/R/differential_expression.R#L663
  
test_data.use <- pseudo_sobj@assays$RNA$counts %>%
  as.data.frame() %>%
  select(contains("CD14 Mono"))

test_group.info <- data.frame(sample = colnames(test_data.use)) %>%
  mutate(group = str_extract(string = sample, pattern = "CTRL|STIM")) %>%
  column_to_rownames(var = "sample")
  
dds1 <- DESeq2::DESeqDataSetFromMatrix(
  countData = test_data.use,
  colData = test_group.info,
  design = ~ group
  )

dds1 <- DESeq2::estimateSizeFactors(object = dds1)
dds1 <- DESeq2::estimateDispersions(object = dds1, fitType = "local")
dds1 <- DESeq2::nbinomWaldTest(object = dds1)

res <- DESeq2::results(
  object = dds1,
  contrast = c("group", "STIM", "CTRL"),
  alpha = 0.05
)

res %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  filter(gene %in% c("CFB", "USP18", "ABTB2"))

de_test1 <- FindMarkers(object = pseudo_sobj,
                        ident.1 = "CD14 Mono_STIM",
                        ident.2 = "CD14 Mono_CTRL",
                        test.use = "DESeq2") %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  filter(gene %in% c("CFB", "USP18", "ABTB2"))

# notice that as mentioned the pvalue si the metric to compare.

# replicate the log2FC ----------------------------------------------------
mono.de <- FindMarkers(sobj, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", verbose = FALSE)
bulk.mono.de <- FindMarkers(object = pseudo_sobj, 
                            ident.1 = "CD14 Mono_STIM", 
                            ident.2 = "CD14 Mono_CTRL",
                            test.use = "DESeq2")

df_full <- mono.de %>%
  rownames_to_column("gene") %>%
  left_join(bulk.mono.de %>%
              rownames_to_column("gene"),
            by = c("gene"),
            suffix = c(".sc",".pBulk"))

# try to understand why the estimates of the logFC are more conservative in the pseudobulk analysis.
df_full %>%
  filter(avg_log2FC.sc>5,avg_log2FC.pBulk<1) %>%
  select(gene,avg_log2FC.sc,avg_log2FC.pBulk,pct.1.sc,pct.2.sc) %>%
  head()

df_full %>%
  filter(abs(avg_log2FC.sc - avg_log2FC.pBulk) < 0.1) %>%
  select(gene,avg_log2FC.sc,avg_log2FC.pBulk,pct.1.sc,pct.2.sc) %>%
  head()

# pbulk
# Get the normalized expression data (log-normalized counts)
expr_matrix <- pseudo_sobj@assays$RNA$data  # This is log-normalized data (ln counts)

# Get the cell names for each group
Idents(pseudo_sobj) <- "celltype.stim"

cells_1 <- WhichCells(pseudo_sobj, idents = "CD14 Mono_STIM")
cells_2 <- WhichCells(pseudo_sobj, idents = "CD14 Mono_CTRL")

# cells_1 <- colnames(expr_matrix) %>%
#   str_subset(pattern = "CD14 Mono") %>%
#   str_subset(pattern = "STIM")
# 
# cells_2 <- colnames(expr_matrix) %>%
#   str_subset(pattern = "CD14 Mono") %>%
#   str_subset(pattern = "CTRL")

# x <- expr_matrix[, cells_1]
avg_expr_1 <- log((rowSums(expm1(expr_matrix[, cells_1])) + 1)/NCOL(expr_matrix[, cells_1]), base = 2)

# x <- expr_matrix[, cells_2]
avg_expr_2 <- log((rowSums(expm1(expr_matrix[, cells_2])) + 1)/NCOL(expr_matrix[, cells_2]), base = 2)

# Compute the natural log fold change (lnFC)
# log2FC <- log(avg_expr_1 + 1) - log(avg_expr_2 + 1)  # Adding 1 avoids log(0)
# log2FC <- avg_expr_1 - avg_expr_2

data.frame(avg_expr_1,avg_expr_2) %>%
  rownames_to_column("gene") %>%
  mutate(log2FC = avg_expr_1 - avg_expr_2) %>%
  filter(gene %in% c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1"))

# compare
df_full %>%
  filter(gene %in% c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1"))

data.frame(s1_STIM = rowSums(expm1(expr_matrix[, cells_1]))[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s1_CTRL = rowSums(expm1(expr_matrix[, cells_2]))[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s2_STIM = (rowSums(expm1(expr_matrix[, cells_1])) + 1)[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s2_CTRL = (rowSums(expm1(expr_matrix[, cells_2])) + 1)[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s3_STIM = ((rowSums(expm1(expr_matrix[, cells_1])) + 1)/NCOL(expr_matrix[, cells_1]))[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s3_CTRL = ((rowSums(expm1(expr_matrix[, cells_2])) + 1)/NCOL(expr_matrix[, cells_2]))[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s4_STIM = log((rowSums(expm1(expr_matrix[, cells_1])) + 1)/NCOL(expr_matrix[, cells_1]),base=2)[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s4_CTRL = log((rowSums(expm1(expr_matrix[, cells_2])) + 1)/NCOL(expr_matrix[, cells_2]),base=2)[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")]) %>%
  mutate(log2FC_s1 = log2(s1_STIM/s1_CTRL),
         log2FC_s2 = log2(s2_STIM/s2_CTRL),
         log2FC_s3 = log2(s3_STIM/s3_CTRL),
         log2FC = s4_STIM - s4_CTRL)

# sc
# Get the normalized expression data (log-normalized counts)
expr_matrix2 <- sobj@assays$RNA$data  # This is log-normalized data (ln counts)

# Get the cell names for each group
Idents(sobj) <- "celltype.stim"

cells_12 <- WhichCells(sobj, idents = "CD14 Mono_STIM")
cells_22 <- WhichCells(sobj, idents = "CD14 Mono_CTRL")

# x <- expr_matrix[, cells_1]
avg_expr_12 <- log((rowSums(expm1(expr_matrix2[, cells_12])) + 1)/NCOL(expr_matrix2[, cells_12]), base = 2)

# x <- expr_matrix[, cells_2]
avg_expr_22 <- log((rowSums(expm1(expr_matrix2[, cells_22])) + 1)/NCOL(expr_matrix2[, cells_22]), base = 2)

# Compute the natural log fold change (lnFC)
# log2FC <- log(avg_expr_1 + 1) - log(avg_expr_2 + 1)  # Adding 1 avoids log(0)
# log2FC2 <- avg_expr_12 - avg_expr_22

data.frame(avg_expr_12,avg_expr_22) %>%
  rownames_to_column("gene") %>%
  mutate(log2FC = avg_expr_12 - avg_expr_22) %>%
  filter(gene %in% c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1"))

df_full %>%
  filter(gene %in% c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1"))

data.frame(s1_STIM = rowSums(expm1(expr_matrix2[, cells_12]))[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s1_CTRL = rowSums(expm1(expr_matrix2[, cells_22]))[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s2_STIM = (rowSums(expm1(expr_matrix2[, cells_12])) + 1)[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s2_CTRL = (rowSums(expm1(expr_matrix2[, cells_22])) + 1)[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s3_STIM = ((rowSums(expm1(expr_matrix2[, cells_12])) + 1)/NCOL(expr_matrix2[, cells_12]))[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s3_CTRL = ((rowSums(expm1(expr_matrix2[, cells_22])) + 1)/NCOL(expr_matrix2[, cells_22]))[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s4_STIM = log((rowSums(expm1(expr_matrix2[, cells_12])) + 1)/NCOL(expr_matrix2[, cells_12]),base=2)[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")],
           s4_CTRL = log((rowSums(expm1(expr_matrix2[, cells_22])) + 1)/NCOL(expr_matrix2[, cells_22]),base=2)[c("CFB", "USP18", "ABTB2","OAS1","ISG20","MX1")]) %>%
  mutate(log2FC_s1 = log2(s1_STIM/s1_CTRL),
         log2FC_s2 = log2(s2_STIM/s2_CTRL),
         log2FC_s3 = log2(s3_STIM/s3_CTRL),
         log2FC = s4_STIM - s4_CTRL)

# simulate a dataset ------------------------------------------------------
library(Seurat)
library(Matrix)

# Set seed for reproducibility
set.seed(42)

# Simulate a small gene expression matrix (raw counts)
genes <- paste0("Gene", 1:10)  # 10 genes
cells <- paste0("Cell", 1:2000)

# Create a raw count matrix with Poisson-distributed values
expr_counts <- matrix(rpois(2000*10, lambda = 10),
                      nrow = 10,
                      ncol = 2000,
                      dimnames = list(genes, cells))

# simulate the addition of a small subset of cells with very low expression of genes

# Convert to a Seurat object
seurat_obj <- CreateSeuratObject(counts = expr_counts)

# create some metadata
seurat_obj$stim <- rep(c("STIM", "CTRL"), each = 1000)
seurat_obj$sample <- rep(paste0("sample",1:10))

# Simulate a gene with higher expression in "STIM" group
gene_of_interest <- "Gene1"

# Increase expression of "Gene1" in "STIM" cells
expr_counts[gene_of_interest, seurat_obj$stim == "STIM"] <- 
  rpois(sum(seurat_obj$stim == "STIM"), lambda = 50)  # Increase lambda for STIM

# Keep expression lower for "CTRL" cells
expr_counts[gene_of_interest, seurat_obj$stim == "CTRL"] <- 
  rpois(sum(seurat_obj$stim == "CTRL"), lambda = 10)  # Maintain lower expression in CTRL

# Update the Seurat object with the modified expression matrix
seurat_obj@assays$RNA$counts <- expr_counts

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)


#####
table(seurat_obj$sample,seurat_obj$stim)

# 
test_expr_matrix <- seurat_obj@assays$RNA$data  # This is log-normalized data (ln counts)

# Get the cell names for each group
# Get the cell names for each group
Idents(seurat_obj) <- "stim"

test_cells_1 <- WhichCells(seurat_obj, idents = "STIM")
test_cells_2 <- WhichCells(seurat_obj, idents = "CTRL")

# x <- expr_matrix[, cells_1]
test_avg_expr_1 <- log((rowSums(expm1(test_expr_matrix[, test_cells_1])) + 1)/NCOL(test_expr_matrix[, test_cells_1]), base = 2)

# x <- expr_matrix[, cells_2]
test_avg_expr_2 <- log((rowSums(expm1(test_expr_matrix[, test_cells_2])) + 1)/NCOL(test_expr_matrix[, test_cells_2]), base = 2)

# Compute the natural log fold change (lnFC)
# log2FC <- log(avg_expr_1 + 1) - log(avg_expr_2 + 1)  # Adding 1 avoids log(0)
test_log2FC <- data.frame(log2FC = test_avg_expr_1 - test_avg_expr_2)


# simulate the pbulk
seurat_obj_pbulk <- AggregateExpression(seurat_obj, assays = "RNA", return.seurat = T, group.by = c("stim", "sample"))

# 
test_expr_matrix2 <- seurat_obj_pbulk@assays$RNA$data  # This is log-normalized data (ln counts)

# Get the cell names for each group
# Get the cell names for each group
Idents(seurat_obj_pbulk) <- "stim"

test_cells_12 <- WhichCells(seurat_obj_pbulk, idents = "STIM")
test_cells_22 <- WhichCells(seurat_obj_pbulk, idents = "CTRL")

# x <- expr_matrix[, cells_1]
test_avg_expr_12 <- log((rowSums(expm1(test_expr_matrix2[, test_cells_12])) + 1)/NCOL(test_expr_matrix2[, test_cells_12]), base = 2)

# x <- expr_matrix[, cells_2]
test_avg_expr_22 <- log((rowSums(expm1(test_expr_matrix2[, test_cells_22])) + 1)/NCOL(test_expr_matrix2[, test_cells_22]), base = 2)

# Compute the natural log fold change (lnFC)
# log2FC <- log(avg_expr_1 + 1) - log(avg_expr_2 + 1)  # Adding 1 avoids log(0)
test_log2FC2 <- data.frame(log2FC = test_avg_expr_12 - test_avg_expr_22)

test_log2FC2

# -------------------------------------------------------------------------


