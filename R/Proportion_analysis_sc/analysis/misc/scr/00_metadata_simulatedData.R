# AIM ---------------------------------------------------------------------
# sample script to explain the problem with proportion analysis

# libraries ---------------------------------------------------------------
library(tidyverse)

# simulate the data -------------------------------------------------------
set.seed(12345)

# cluster metadata
cluster <- rep(paste0("clu ",c(1:5)),12)

# sample metadata
sample <- rep(LETTERS[1:12],each = 5)

# sample treat
treat <- c(rep("ctrl",30),rep("treat",30))

# build the total dataset
df <- data.frame(count = round(runif(60,min = 500,max = 800),digits = 0)) %>%
  mutate(cluster = cluster) %>%
  mutate(sample = sample) %>%
  mutate(treat = treat)

# simulate an increase in one specific cluster only
df_test <- df %>%
  mutate(count_test = case_when(cluster == "clu 1" & treat == "treat" ~ count + 600,
                                T ~ count))

# wrangling ---------------------------------------------------------------
# calculate the proportion
df_tot <- df_test %>%
  group_by(sample,treat) %>%
  # regular condition
  mutate(total = sum(count),
         prop = count/total) %>%
  # test condition
  mutate(total_test = sum(count_test),
         prop_test = count_test/total_test)
  

# show the proportions
df_tot %>%
  dplyr::select(cluster,sample,treat,prop,prop_test) %>%
  pivot_longer(names_to = "var_name",values_to = "prop",c(prop,prop_test)) %>%
  ggplot(aes(x=cluster,y=prop,col=treat)) +
  geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),shape=1,alpha=0.8)+
  facet_wrap(~var_name,ncol = 1) +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave("../../out/plot/boxplot_proportions.pdf",height = 6,width = 6)  

# save the object ---------------------------------------------------------
df_tot %>%
  write_tsv("../../out/table/sample_metadata_simulatedData.tsv")
