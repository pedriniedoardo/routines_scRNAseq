(DimPlot(data.combined,group.by = "test_cellid",label = T)+ggtitle("so_ref"))+
  (DimPlot(data.combined2,group.by = "test_cellid",label = T)+ggtitle("so_test"))
  

# change color of the plot ------------------------------------------------
#
test1 <- plotNhoodGraphDA(milo, da_results, alpha=0.05)
test2 <- plotNhoodGraphDA(milo2, da_results2, alpha=0.05)

test2+scale_fill_viridis_c(option = "turbo")

test2+scale_fill_gradient2(
  low = "blue",
  mid = "white",
  high = "red",
  midpoint = 0,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill"
)

test2+scale_fill_gradientn(colours = viridis::viridis(20),limits = c(-2,2),oob = scales::squish,name = 'logFC')
library(circlize)
color_values <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))(seq(from=-2,to=2,by=0.1))

(test1+scale_fill_gradientn(limits = c(-2,2),oob = scales::squish,name = 'logFC',colours = color_values)+ggtitle("so_ref"))+
(test2+scale_fill_gradientn(limits = c(-2,2),oob = scales::squish,name = 'logFC',colours = color_values)+ggtitle("so_test"))

test2+scale_fill_gradient2(
  low = "blue",
  mid = "white",
  high = "red",
  midpoint = 0,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill"
)

test2+scale_fill_gradient2(
  low = "blue",
  mid = "white",
  high = "red",
  midpoint = 0,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill"
)

test2$data %>% 
  ggplot(aes(x=x,y=y))+
  geom_point(aes(fill=colour_by,size = size),shape=21,alpha=0.5)+scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  )+theme_bw()
ggsave("out/image/milo2R_test_update2.pdf",width = 7,height = 6)

# Inspecting DA testing results -------------------------------------------
# We can start inspecting the results of our DA analysis from a couple of standard diagnostic plots. We first inspect the distribution of uncorrected P values, to verify that the test was balanced.
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)+theme_bw()+ggtitle("so_ref")+ggplot(da_results2, aes(PValue)) + geom_histogram(bins=50)+theme_bw()+ggtitle("so_test")

p1 <- ggplot(da_results2, aes(logFC, -log10(SpatialFDR))) + 
  geom_point(shape = 1,alpha=0.5) +
  ## Mark significance threshold (5% FDR)
  geom_hline(yintercept = -log10(0.05),linetype="dashed",col="gray") +
  theme_bw()+
  theme(strip.background = element_blank())+
  ggtitle("so_test")

p2 <- ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point(shape = 1,alpha=0.5) +
  ## Mark significance threshold (5% FDR)
  geom_hline(yintercept = -log10(0.05),linetype="dashed",col="gray") +
  theme_bw()+
  theme(strip.background = element_blank())+
  ggtitle("so_ref")

p2+p1

# To visualize DA results relating them to the embedding of single cells, we can build an abstracted graph of neighbourhoods that we can superimpose on the single-cell embedding. Here each node represents a neighbourhood, while edges indicate how many cells two neighbourhoods have in common. Here the layout of nodes is determined by the position of the index cell in the UMAP embedding of all single-cells. The neighbourhoods displaying significant DA are colored by their log-Fold Change.
milo2 <- buildNhoodGraph(milo2)

## Plot single-cell UMAP
umap_pl2 <- plotReducedDim(milo2, dimred = "UMAP", colour_by="test_disease", text_by = "test_cellid", 
                           text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl2 <- plotNhoodGraphDA(milo2, da_results2, layout="UMAP",alpha=0.98) 
umap_pl2 + nh_graph_pl2 + plot_layout(guides="collect")

# We might also be interested in visualizing whether DA is particularly evident in certain cell types. To do this, we assign a cell type label to each neighbourhood by finding the most abundant cell type within cells in each neighbourhood. We can label neighbourhoods in the results data.frame using the function annotateNhoods. This also saves the fraction of cells harbouring the label.
da_results2 <- annotateNhoods(milo2, da_results2, coldata_col = "test_cellid")
head(da_results2)

# Then we visualize the test results with a volcano plot (remember that each point here represents a neighbourhood, not a cell).
p11 <- ggplot(da_results2, aes(logFC, -log10(SpatialFDR))) + 
  geom_point(shape = 1,alpha=0.5) +
  ## Mark significance threshold (5% FDR)
  geom_hline(yintercept = -log10(0.05),linetype="dashed",col="gray") +
  theme_bw()+
  theme(strip.background = element_blank())+
  facet_wrap(~test_cellid)+
  ggtitle("so_test")

p21 <- ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point(shape = 1,alpha=0.5) +
  ## Mark significance threshold (5% FDR)
  geom_hline(yintercept = -log10(0.05),linetype="dashed",col="gray") +
  theme_bw()+
  theme(strip.background = element_blank())+
  facet_wrap(~test_cellid)+
  ggtitle("so_ref")

p21+p11

# While neighbourhoods tend to be homogeneous, we can define a threshold for celltype_fraction to exclude neighbourhoods that are a mix of cell types.
ggplot(da_results2, aes(test_cellid_fraction)) + geom_histogram(bins=50)

da_results2$test_cellid <- ifelse(da_results2$test_cellid_fraction < 0.7, "Mixed", da_results2$test_cellid)
# Now we can visualize the distribution of DA Fold Changes in different cell types
test2 <- plotDAbeeswarm(da_results2, group.by = "test_cellid",alpha = 1)

plotDAbeeswarm

test2$data %>%
  group_by(group_by) %>%
  summarise(med = median(pos_x)) %>%
  arrange(med) %>%
  pull(group_by)

p13 <- test2$data %>%
  mutate(color = case_when(FDR<0.05~"sig",
                           T~"n.s.")) %>%
  arrange(color) %>%
  ggplot(aes(x=logFC,y=group_by,color = color)) +
  ggbeeswarm::geom_quasirandom(alpha=0.5)+
  theme_bw()+
  scale_color_manual(values = c("gray","red"))+
  scale_x_continuous(limits = c(-2.5,2.5))+
  ggtitle("so_test")

p23 <- test$data %>%
  mutate(color = case_when(FDR<0.05~"sig",
                           T~"n.s.")) %>%
  arrange(color) %>%
  ggplot(aes(x=logFC,y=group_by,color = color)) +
  ggbeeswarm::geom_quasirandom(alpha=0.5)+
  theme_bw()+
  scale_color_manual(values = c("gray","red"))+
  scale_x_continuous(limits = c(-2.5,2.5))+
  ggtitle("so_ref")+
  theme(legend.position = "none")

p23+p13

# This is already quite informative: we can see that certain early development cell types, such as epiblast and primitive streak, are enriched in the earliest time stage, while others are enriched later in development, such as ectoderm cells. Interestingly, we also see plenty of DA neighbourhood with a mixed label. This could indicate that transitional states show changes in abundance in time.
