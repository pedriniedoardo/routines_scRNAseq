. /home/edo/micromamba/bin/activate
conda activate env_cellbender_pip

# cellbender -v

cellbender remove-background \
     --input data/tiny_raw_feature_bc_matrix.h5ad \
     --output out/GPU_pip/tiny_output.h5 \
     --expected-cells 500 \
     --total-droplets-included 2000 \
     --cuda