renv::install("SeuratObject")
renv::install("Seurat")
renv::install("bnprks/BPCells")
renv::install("satijalab/seurat-data")
renv::install("R.utils")
renv::install("satijalab/seurat-wrappers")
renv::install("robustbase")

renv::install("bioc::GenomeInfoDb")
# renv::install("bioc::BSgenome.Hsapiens.UCSC.hg38")
# renv::install("bioc::EnsDb.Hsapiens.v86")
# renv::install("bioc::glmGamPoi")
# renv::install("bioc::JASPAR2020")
renv::install("mojaveazure/seurat-disk")
# renv::install("bioc::TFBSTools")

renv::install("bioc::GenomicRanges")
renv::install("bioc::Rsamtools")
renv::install("stuart-lab/signac")
renv::install("satijalab/azimuth")

renv::install("tidyverse")

# -------------------------------------------------------------------------
library(Seurat)
library(BPCells)
library(SeuratData)
library(presto)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(Azimuth)
library(tidyverse)


