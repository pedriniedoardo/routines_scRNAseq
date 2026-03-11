#!/bin/bash
. /home/edo/micromamba/bin/activate
conda activate env_MQuad

# Add --chrom if you only want to genotype specific chromosomes, e.g., 1,2, or chrMT
cellsnp-lite -s /mnt/SSD01/training/routines/scRNAseq/MQuad/data/data_test_cellsnplite/cellSNP_testdata_10x/demux.B.lite.bam \
-b /mnt/SSD01/training/routines/scRNAseq/MQuad/data/data_test_cellsnplite/cellSNP_testdata_10x/demux.B.barcodes.400.tsv \
-O /mnt/SSD01/training/routines/scRNAseq/MQuad/test/out/example_cellsnplite_02_mito \
-p 10 \
--minMAF 0 \
--minCOUNT 10 \
--chrom MT \
--gzip