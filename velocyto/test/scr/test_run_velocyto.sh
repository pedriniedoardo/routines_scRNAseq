. /home/edo/micromamba/bin/activate
conda activate env_velocyto

# the output will be created inside the cellranger output folder as velocyto/name_sample.loom

velocyto run10x \
-m /mnt/SSD01/training/routines/scRNAseq/velocyto/ref/repeat_mask_GRCh38.gtf \
/mnt/SSD01/training/routines/scRNAseq/velocyto/test/out_cellranger/cellranger901/test_pbmc \
/home/edo/Documents/reference/cellranger/reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf