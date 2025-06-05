library(Seurat)
library(ggplot2)
library(patchwork)

library(SeuratData)
library(Azimuth)


reference <- readRDS("../data/pbmc_multimodal_2023.rds")

set.seed(2144)
ref_meta <- reference@meta.data %>%
  rownames_to_column("barcode") %>%
  sample_n(5000)

id_test <- reference@meta.data %>%
  rownames_to_column("barcode") %>%
  mutate(test = barcode %in% ref_meta$barcode) %>%
  pull(test)

# update the meta
reference$test <- id_test

reference_subset <- subset(reference,subset = test == 1)

saveRDS(reference_subset,"../out/object/pbmc_multimodal_2023_subset.rds")

AvailableData()
InstallData('pbmc3k')

pbmc <- LoadData('pbmc3k')

LoadData('pbmcsca')
pbmc3k <- UpdateSeuratObject(pbmc3k)

# InstallData("pbmcMultiome")
# InstallData("pbmcMultiome.SeuratData")

LoadData("pbmcMultiome")

# download.file(url = "http://seurat.nygenome.org/src/contrib/pbmcMultiome.SeuratData_0.1.4.tar.gz",destfile = "../data/pbmcMultiome.SeuratData_0.1.4.tar.gz")
download.file(url = "http://seurat.nygenome.org/src/contrib/pbmcMultiome.SeuratData_0.1.2.tar.gz",destfile = "../data/pbmcMultiome.SeuratData_0.1.2.tar.gz")
install.packages("../data/pbmcMultiome.SeuratData_0.1.2.tar.gz")

LoadData("pbmcMultiome.SeuratData")


LoadFileInput("renv/library/R-4.3/x86_64-pc-linux-gnu/pbmc3k.SeuratData/data/")
pbmc <- LoadData('pbmc3k')
reference <- LoadReference(path = "renv/library/R-4.3/x86_64-pc-linux-gnu/pbmcMultiome.SeuratData/data")
reference <- LoadData("renv/library/R-4.3/x86_64-pc-linux-gnu/pbmc3k.SeuratData/")

Azimuth::LoadReference(path = "renv/library/R-4.3/x86_64-pc-linux-gnu/pbmcMultiome.SeuratData/data", "default")

reference <- LoadData("renv/library/R-4.3/x86_64-pc-linux-gnu/pbmcMultiome.SeuratData/data")#Could not find dataset 'pbmcMultiome'
data(package = "pbmcMultiome.SeuratData")
