. /home/edo/micromamba/bin/activate;
conda activate env_cellranger901

cellranger count --id=test_neuron_auto \
        --transcriptome=/home/edo/Documents/reference/cellranger/reference/refdata-gex-GRCm39-2024-A \
        --fastqs=/mnt/SSD01/training/routines/scRNAseq/cellranger/test_annotation/data/neuron_1k_v3_fastqs \
        --sample=neuron_1k_v3 \
        --localcores=8 \
        --create-bam=false \
        --cell-annotation-model=auto \
        --localmem=64 \
        --output-dir=/mnt/SSD01/training/routines/scRNAseq/cellranger/test_annotation/out/cellranger901/test_neuron_auto