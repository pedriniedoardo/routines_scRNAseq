# AIM ---------------------------------------------------------------------
# reccommended processing for the graph object in case of big rds objects.
# see the discussion of github for more details:
# https://github.com/MarioniLab/miloR/issues/351#issuecomment-2462720639

# libraries ---------------------------------------------------------------
library(miloR)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)

# read in data ------------------------------------------------------------
# read in a sample big dataset
ref <- readRDS(file = "../../data/harmonyHO.rds")

# confirm the detaul slot for the data is RNA
DefaultAssay(ref)

# wrangling ---------------------------------------------------------------
# extract the graph object
graph_obj <- ref@graphs$RNA_snn

# check the size and the stricture of the object
graph_obj.size <- object.size(graph_obj)
print(graph_obj.size,units = "Mb")
# confirm it is a sparse matrix
str(graph_obj)

# make it as a sparse matrix
graph_obj_fix <- as(graph_obj, "dgCMatrix")

min(graph_obj_fix)
max(graph_obj_fix)
sum(is.na(graph_obj_fix))

sum(graph_obj_fix)

# Follow the reccomandation suggested on github and make it as a binary matrix
graph_obj_fix2 <- graph_obj_fix
graph_obj_fix2@x <- rep(x = 1,length(graph_obj_fix2@x))
graph_obj_fix.size <- object.size(graph_obj_fix2)
print(graph_obj_fix.size,units = "Mb")

# build the adiacency matrix using binary as T. the RAM should not spike anymore
test <- buildFromAdjacency(graph_obj_fix2, k=10,is.binary = T)
test.size <- object.size(test)
print(test.size,units = "Mb")

# create the milo object --------------------------------------------------
# create the milo object using the graph object calculated in seurat

# build the milo object
sce3 <- as.SingleCellExperiment(ref)
milo3 <- Milo(sce3)

# add the graph generated to the object
miloR::graph(milo3) <- miloR::graph(test)

# milo processing ---------------------------------------------------------
# confirm the simentional reduction is correctly interpreted
plotUMAP(milo3)

# confirm the metadata are loaded
colData(milo3)

# confirm meta and data are correctly loaded
plotUMAP(milo3,colour_by="ident")

# confirm the community are generated using the graph object
milo3 <- makeNhoods(milo3, prop = 0.1, k = 10, d=30, refined = TRUE)
plotNhoodSizeHist(milo3)
