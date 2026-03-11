#!/bin/bash
. /home/edo/micromamba/bin/activate
conda activate env_MQuad

mquad --vcfData /mnt/SSD01/training/routines/scRNAseq/MQuad/data/example.vcf.gz \
-o /mnt/SSD01/training/routines/scRNAseq/MQuad/test/out/example_mquad_01 \
-p 5