# AIM ---------------------------------------------------------------------
# Test the scProportionTest vignette

# libraries ---------------------------------------------------------------
library("scProportionTest")
library("tidyverse")

# read in the data --------------------------------------------------------
seurat_data <- system.file("extdata", "example_data.RDS", package = "scProportionTest")
seurat_data <- readRDS(seurat_data)

# update the object
seurat_data_update <- Seurat::UpdateSeuratObject(seurat_data)

# sample processing -------------------------------------------------------
# the function extract the meta from the object
prop_test <- sc_utils(seurat_data_update)

# characterize the samples
prop_test@meta_data %>%
  group_by(custom_clusters) %>%
  summarise(n = n())

prop_test@meta_data %>%
  group_by(orig.ident) %>%
  summarise(n = n())

# the function run the test camparing
# arguments
# sc_utils_obj: sc_utils object
# cluster_identity: Column that has cluster names
# sample_1: First sample to compare (ie. control)
# sample_2: Sample to compare to first sample (ie. treatment)
# sample_identity: Column that has sample names
# n_permutations: Number of permutations
prop_test_results <- permutation_test(prop_test,
                                      cluster_identity = "custom_clusters",
                                      sample_1 = "HT29_EV",
                                      sample_2 = "HT29_LSD1_KD",
                                      sample_identity = "orig.ident")

# plot the data
permutation_plot(prop_test_results)
