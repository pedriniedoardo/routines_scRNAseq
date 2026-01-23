. /home/edo/micromamba/bin/activate;
conda activate env_cellranger901

cellranger count --id=test_pbmc \
        --transcriptome=/home/edo/Documents/reference/cellranger/reference/refdata-gex-GRCh38-2024-A \
        --fastqs=/mnt/SSD01/training/routines/scRNAseq/velocyto/ref/fastq/connect_5k_pbmc_NGSC3_ch1_gex_2_tiny \
        --sample=connect_5k_pbmc_NGSC3_ch1_gex_2 \
        --localcores=8 \
        --create-bam=true \
        --cell-annotation-model=auto \
        --localmem=64 \
        --output-dir=/mnt/SSD01/training/routines/scRNAseq/velocyto/test/out_cellranger/cellranger901/test_pbmc