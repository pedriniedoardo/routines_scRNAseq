#!/bin/bash
. /home/edo/micromamba/bin/activate
conda activate env_MQuad

cellsnp-lite -s /mnt/SSD01/training/routines/scRNAseq/MQuad/data/data_test_cellsnplite/cellSNP_testdata_10x/demux.B.lite.bam \
-b /mnt/SSD01/training/routines/scRNAseq/MQuad/data/data_test_cellsnplite/cellSNP_testdata_10x/demux.B.barcodes.400.tsv \
-O /mnt/SSD01/training/routines/scRNAseq/MQuad/test/out/example_cellsnplite_01 \
-R /mnt/SSD01/training/routines/scRNAseq/MQuad/data/data_test_cellsnplite/cellSNP_testdata_10x/genome1K.subset.hg19.vcf.gz \
-p 10 \
--minMAF 0.1 \
--minCOUNT 20 \
--gzip