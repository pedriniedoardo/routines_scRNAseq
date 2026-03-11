#!/bin/bash
. /home/edo/micromamba/bin/activate
conda activate env_MQuad

vireo -c /mnt/SSD01/training/routines/scRNAseq/MQuad/data/data_test_vireo/cellSNP_mat \
-o /mnt/SSD01/training/routines/scRNAseq/MQuad/test/out/example_vireo_mode_01 \
-N 4 \
--randSeed 2