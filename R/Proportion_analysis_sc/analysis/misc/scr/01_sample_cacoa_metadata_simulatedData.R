# AIM ---------------------------------------------------------------------
# run a statistical test on the proportion differences across clusters, between conditions

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(cacoa)
library(conos)
library(sccore)
library(coda.base)
library(psych)
library(quadprog)
library(patchwork)

# read in the data --------------------------------------------------------
# read in the summary simualted metadata
df <- read_tsv("../../out/table/sample_metadata_simulatedData.tsv")

# generate a fake metadata based on the summary count
# ref data
meta_ref <- df %>%
  mutate(col_split = paste(sample,cluster,treat,sep = "_")) %>%
  split(f = .$col_split) %>%
  lapply(function(x){
    # determine the number of replications of the entry in the metada
    n_row <- x$count
    
    # define the columns
    cluster <- rep(x$cluster,n_row)
    sample <- rep(x$sample,n_row)
    treat <- rep(x$treat,n_row)
    
    # build the fake metadata
    meta <- data.frame(sample,treat,cluster)
    return(meta)
  }) %>%
  bind_rows() %>%
  rownames_to_column("rowname") %>%
  mutate(rowname2 = paste0("Cell_",rowname)) %>%
  column_to_rownames("rowname2")

# test stimulation
meta_test <- df %>%
  mutate(col_split = paste(sample,cluster,treat,sep = "_")) %>%
  split(f = .$col_split) %>%
  lapply(function(x){
    # determine the number of replications of the entry in the metada
    n_row <- x$count_test
    
    # define the columns
    cluster <- rep(x$cluster,n_row)
    sample <- rep(x$sample,n_row)
    treat <- rep(x$treat,n_row)
    
    # build the fake metadata
    meta <- data.frame(sample,treat,cluster)
    return(meta)
  }) %>%
  bind_rows() %>%
  rownames_to_column("rowname") %>%
  mutate(rowname2 = paste0("Cell_",rowname)) %>%
  column_to_rownames("rowname2")

# check dimensions --------------------------------------------------------
# confirm the numbers from the reference dataset
head(meta_ref)
dim(meta_ref)
sum(df$count)

head(meta_test)
dim(meta_test)
sum(df$count_test)

# wrangling ---------------------------------------------------------------
# create a fake seurat object with the current metadate
mat_ref <- as.sparse(matrix(0,nrow = 2,ncol = dim(meta_ref)[1]))
obj_ref <- CreateSeuratObject(counts = mat_ref,meta.data = meta_ref)

mat_test <- as.sparse(matrix(0,nrow = 2,ncol = dim(meta_test)[1]))
obj_test <- CreateSeuratObject(counts = mat_test,meta.data = meta_test)

# run cacoa ---------------------------------------------------------------
# sample.groups: vector with condition labels per sample named with sample ids
sample.groups_ref <- meta_ref$treat
names(sample.groups_ref) <- meta_ref$sample

sample.groups_test <- meta_test$treat
names(sample.groups_test) <- meta_test$sample

# cell.groups: cell type annotation vector named by cell ids
cell.groups_ref <- meta_ref$cluster
names(cell.groups_ref) <- rownames(meta_ref)

cell.groups_test <- meta_test$cluster
names(cell.groups_test) <- rownames(meta_test)

# sample.per.cell: vector with sample labels per cell named with cell ids
sample.per.cell_ref <- meta_ref$sample
names(sample.per.cell_ref) <- rownames(meta_ref)

sample.per.cell_test <- meta_test$sample
names(sample.per.cell_test) <- rownames(meta_test)

# ref.level: id of the condition, corresponding to the reference (i.e. control)
table(sample.groups_ref)
ref.level <- "ctrl"

# target.level: id of the condition, corresponding to the target (i.e. case)
target.level <- "treat"

# creat the cacoa object
# cao_ref <- Cacoa$new(data.object = obj_ref,
#                      sample.groups=sample.groups_ref,
#                      cell.groups=cell.groups_ref,
#                      sample.per.cell=sample.per.cell_ref,
#                      ref.level=ref.level,
#                      target.level=target.level)
# cao_test <- Cacoa$new(data.object = obj_test,
#                       sample.groups=sample.groups_test,
#                       cell.groups=cell.groups_test,
#                       sample.per.cell=sample.per.cell_test,
#                       ref.level=ref.level,
#                       target.level=target.level)

# No expression data
# Cacoa can be ran without any expression data by passing NULL instead of a data object:
# In this case, only compositional analyses will be available.
cao_ref <- Cacoa$new(data.object = NULL,
                     sample.groups=sample.groups_ref,
                     cell.groups=cell.groups_ref,
                     sample.per.cell=sample.per.cell_ref,
                     ref.level=ref.level,
                     target.level=target.level)

cao_test <- Cacoa$new(data.object = NULL,
                      sample.groups=sample.groups_test,
                      cell.groups=cell.groups_test,
                      sample.per.cell=sample.per.cell_test,
                      ref.level=ref.level,
                      target.level=target.level)

# run cacoa Cluster-based changes -----------------------------------------
# The fastest way to visualize changes in the dataset is to show, what cell types changed their abundance and expression patterns.

# Estimate cluster-based changes
cao_ref$estimateCellLoadings()
cao_test$estimateCellLoadings()

# Plot compositional changes
cao_ref$plotCellLoadings(show.pvals=T)
cao_test$plotCellLoadings(show.pvals=T)

# order_cluster <- attr(cao_ref$test.results$coda$loadings, "dimnames")[[1]]
# 
# cao_ref$test.results$coda$loadings %>%
#   data.frame() %>%
#   rownames_to_column("cluster") %>%
#   pivot_longer(names_to = "cell",values_to = "loading",-cluster) %>%
#   mutate(cluster = fct_relevel(cluster,rev(order_cluster))) %>%
#   ggplot(aes(x=loading,y=cluster))+geom_boxplot(outlier.shape = NA)+
#   geom_vline(xintercept = 0)

list_plot <- lapply(list(cao_ref=cao_ref,cao_test=cao_test),function(x){
  
  alpha <- 0.01
  palette <- x$cell.groups.palette
  font.size=NULL
  name='coda'
  ordering='pvalue'
  show.pvals=TRUE
  
  loadings <- x$test.results$coda$loadings
  pval <- x$test.results$coda$padj
  ref.level
  target.level
  signif.threshold=0.05
  jitter.alpha=0.1
  # palette=cao_ref$cell.groups.palette
  # show.pvals=T
  plot.theme=theme_bw()
  jitter.size=1
  # ordering=c("pvalue", "loading")
  # ordering=c("pvalue")
  ref.load.level=0
  annotation.x=NULL
  annotation.y=1
  
  
  # replicate plotting ------------------------------------------------------
  # plotCellLoadings <- function(loadings, pval, ref.level, target.level, signif.threshold=0.05, jitter.alpha=0.1,
  #                              palette=NULL, show.pvals=FALSE, plot.theme=theme_get(), jitter.size=1,
  #                              ordering=c("pvalue", "loading"), ref.load.level=0, annotation.x=NULL, annotation.y=1) {
  ordering <- "pvalue"
  xintercept <- ref.load.level
  
  loading.order <- order(abs(rowMeans(loadings)))
  loadings <- loadings[loading.order, ]
  
  if (!is.null(pval)) {
    # if some p-values are the same - then order by mean, therefore a prior sorting is required
    pval <- pval[loading.order]
    
    if (ordering == "pvalue") {
      # additional ordering by p-value
      loadings <- loadings[order(-pval), ]
      pval <- pval[order(-pval)]
    }
    
    # Get significant cells
    n.significant.cells <- sum(pval <= signif.threshold)
  } else {
    n.significant.cells <- 0
  }
  
  # Normalization of loadings
  loadings <- loadings - xintercept
  
  res.ordered <- t(loadings) %>% as.data.frame()
  
  if (is.null(annotation.x)) {
    annotation.x <- max(loadings)
  }
  
  p <- ggplot(stack(res.ordered), aes(y=ind, x=values, fill=factor(ind))) +
    geom_boxplot(notch=TRUE, outlier.shape = NA)
  
  if ((jitter.alpha > 1e-5) && (jitter.size > 1e-5)) {
    p <- p + geom_jitter(alpha=jitter.alpha, size=jitter.size)
  }
  
  p <- p + geom_vline(xintercept=0, color = "grey37") +
    labs(y='', x='loadings') + plot.theme + theme(legend.position = "none") +
    scale_y_discrete(position = "right") + xlim(-1, 1)
  
  # Add text
  p <- p +
    annotate('text', x=-annotation.x, y=annotation.y, label=paste('\u2190', ref.level), hjust='left') +
    annotate('text', x=annotation.x, y=annotation.y, label=paste(target.level, '\u2192'), hjust='right')
  
  if ((n.significant.cells > 0) && (ordering == "pvalue")) {
    p <- p + geom_hline(yintercept=nrow(loadings) - n.significant.cells + 0.5, color='red')
  }
  
  if (!is.null(palette)) p <- p + scale_fill_manual(values=palette)
  
  d <- data.frame(y=-log(pval, base=10), x=names(pval), row=1:length(pval))
  
  p.pval <- ggplot(d, aes(y=reorder(x, row), x=y, fill=factor(x))) +
    geom_bar(stat="identity") +
    geom_vline(xintercept=-log(signif.threshold, base=10)) +
    scale_x_continuous(expand=c(0, 0, 0.02, 0),limits = c(0,2.5)) +
    labs(y='', x='-log10(adj. p-value)') +
    plot.theme +
    theme(legend.position="none") +
    theme(axis.text.y=element_blank())
  
  if (!is.null(palette)) p.pval <- p.pval + scale_fill_manual(values=palette)
  
  p.combo <- cowplot::plot_grid(p, p.pval, nrow=1, rel_widths=c(2, 1))
  return(p.combo)
  # return(list(loading = p,
  #             pvalue = p.pval))
})


list_plot$cao_ref+plot_annotation("meta_ref")
list_plot$cao_test+plot_annotation("meta_test")
