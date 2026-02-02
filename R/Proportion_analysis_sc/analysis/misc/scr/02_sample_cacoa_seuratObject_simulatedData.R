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

# read in the data --------------------------------------------------------
# create a fake seurat object with the current metadate
obj_ref <- readRDS("../../out/object/sobj_total_full_h.rds")
obj_test <- readRDS("../../out/object/sobj_total_remove_h.rds")

meta_ref <- obj_ref@meta.data
meta_test <- obj_test@meta.data

# run cacoa ---------------------------------------------------------------
# sample.groups: vector with condition labels per sample named with sample ids
sample.groups_ref <- meta_ref$test_disease
names(sample.groups_ref) <- meta_ref$test_donor

sample.groups_test <- meta_test$test_disease
names(sample.groups_test) <- meta_test$test_donor

# cell.groups: cell type annotation vector named by cell ids
cell.groups_ref <- meta_ref$test_cellid
names(cell.groups_ref) <- rownames(meta_ref)

cell.groups_test <- meta_test$test_cellid
names(cell.groups_test) <- rownames(meta_test)

# sample.per.cell: vector with sample labels per cell named with cell ids
sample.per.cell_ref <- meta_ref$test_donor
names(sample.per.cell_ref) <- rownames(meta_ref)

sample.per.cell_test <- meta_test$test_donor
names(sample.per.cell_test) <- rownames(meta_test)

# ref.level: id of the condition, corresponding to the reference (i.e. control)
table(sample.groups_ref)
ref.level <- "ctrl"

# target.level: id of the condition, corresponding to the target (i.e. case)
target.level <- "dis"

# creat the cacoa object
cao_ref <- Cacoa$new(data.object = obj_ref,
                     sample.groups=sample.groups_ref,
                     cell.groups=cell.groups_ref,
                     sample.per.cell=sample.per.cell_ref,
                     ref.level=ref.level,
                     target.level=target.level,
                     data.slot = "data")

cao_test <- Cacoa$new(data.object = obj_test,
                      sample.groups=sample.groups_test,
                      cell.groups=cell.groups_test,
                      sample.per.cell=sample.per.cell_test,
                      ref.level=ref.level,
                      target.level=target.level,
                      data.slot = "data")

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
  
  alpha <- 0.001
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
  jitter.alpha=0.01
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


list_plot

list_plot$cao_ref+plot_annotation("cao_ref")
list_plot$cao_test+plot_annotation("cao_test")
