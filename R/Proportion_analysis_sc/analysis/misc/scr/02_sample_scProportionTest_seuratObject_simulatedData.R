# AIM ---------------------------------------------------------------------
# test sc proportion on the simulated metadata

# libraries ---------------------------------------------------------------
library("scProportionTest")
library("tidyverse")
library("scales")

# read in the data --------------------------------------------------------
# create a fake seurat object with the current metadate
obj_ref <- readRDS("../../out/object/sobj_total_full_h.rds")
obj_test <- readRDS("../../out/object/sobj_total_remove_h.rds")

# sample processing -------------------------------------------------------
# the function extract the meta from the object
prop_ref <- sc_utils(obj_ref)
prop_test <- sc_utils(obj_test)

# sanity check
head(prop_ref@meta_data)
head(prop_test@meta_data)

table(prop_ref@meta_data$test_cellid)
table(prop_ref@meta_data$test_disease)
table(prop_ref@meta_data$test_donor)

# the function run the test camparing
# arguments
# sc_utils_obj: sc_utils object
# cluster_identity: Column that has cluster names
# sample_1: First sample to compare (ie. control)
# sample_2: Sample to compare to first sample (ie. treatment)
# sample_identity: Column that has sample names
# n_permutations: Number of permutations
prop_ref_results <- permutation_test(prop_ref,
                                     cluster_identity = "test_cellid",
                                     sample_1 = "ctrl",
                                     sample_2 = "dis",
                                     sample_identity = "test_disease")

prop_test_results <- permutation_test(prop_test,
                                      cluster_identity = "test_cellid",
                                      sample_1 = "ctrl",
                                      sample_2 = "dis",
                                      sample_identity = "test_disease")

# plot the data
# permutation_plot(prop_ref_results,log2FD_threshold = 0.2)
# permutation_plot(prop_test_results,log2FD_threshold = 0.2)

# make a unique plot with both data
df_plot <- bind_rows (res_ref <- prop_ref_results@results$permutation %>%
                        data.frame() %>%
                        mutate(dataset = "prop_ref"),
                      res_test <- prop_test_results@results$permutation %>%
                        data.frame() %>%
                        mutate(dataset = "prop_test")) %>%
  mutate(sig = case_when(FDR < 0.05 & abs(boot_mean_log2FD) > 0.1 ~ "FDR < 0.05 & abs(boot_mean_log2FD) > 0.1",
                         T ~ "n.s"))

# plot
show_col(hue_pal()(1))
df_plot %>%
  mutate(id = paste(dataset,clusters,sep = "|")) %>%
  pivot_longer(names_to = "side",values_to = "x",c(boot_CI_2.5,boot_CI_97.5)) %>%
  ggplot(aes(x=x,y=clusters,group = clusters,col=sig))+
  geom_line()+
  geom_point(data = df_plot,aes(x=boot_mean_log2FD,y=clusters,col=sig)) +
  facet_wrap(~dataset,ncol = 1)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = c(-0.1,+0.1),linetype="dashed")+
  scale_color_manual(values = c(hue_pal()(1),"gray"))+
  xlab("boot_mean_log2FD")
