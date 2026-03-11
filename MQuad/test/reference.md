# set-up
create the environment with all the tools needed
- cellsnp-lite to pileup mtDNA variants from raw .bam file(s).
- MQuad to differentiate informative mtDNA variants from noisy backbground.
- vireoSNP to assign cells to clones based on mtDNA variant profile.
use as reference the yaml file `env_MQuad`

# pile up variants
The instruction below come from the `preprocessing_smd.sh` from the mquad github repo

## SMART-SEQ2
For sequencing data with individual .bam files for each cell + no UMItags/barcodes run `cellsnp-lite` mode2a on bam list. change `--chrom=` to whatever reference genome you aligned to. in this case we use `hg19`.

```
ls *.bam > bam.lst
cellsnp-lite -S bam.lst -i bam.lst -o cellsnp --UMItag None --genotype --gzip --chrom=chrM -p 10
```

## 10X/UMI-BASED
For 10x and UMI-based sequencing data where there is only 1 big .bam + barcodes run `cellsnp-lite` on the `.bam` directly. Change `--chrom=` to whatever reference genome you aligned to, in most cases 10x data are aligned to GrCh38 so the `chr` name is `MT`.

```
cellsnp-lite -s possorted_genome.bam -b barcodes.tsv -o cellsnp --chrom=MT --UMItag Auto --minMAF 0 --minCOUNT 0 --genotype --gzip -p 10
```

The above steps should generate a cell x snp .vcf file (`cellSNP.cells.vcf.gz`), or AD/DP sparse matrices if you did not use the `--genotype` option. Then you should be able to run MQuad on the vcf file.

# run MQuad
MQuad comes with an example dataset for you to test things out. The mtDNA mutations of this dataset are extracted from Ludwig et al, Cell, 2019. It contains 500 background variants, along with 9 variants used in Supp Fig. 2F (and main Fig. 2F). There is also 1 additional variant that is informative but not mentioned in the paper. In total, there are 510 variants in the example dataset.

Run the following command line:

```
mquad --vcfData example/example.vcf.gz -o example_test -p 5
```

or using batch mode tailored for mixture-binomial modelling:

```
mquad --vcfData example/example.vcf.gz -o example_test -p 5 --batchFit 1 --batchSize 5
```

# run vireo
It is generally difficult to identify the number of clones, which is a balance between subclone resolution and analysis reliability. More clones maybe preferred, but there could be higher risk that the subclones are not genuine but rather technical noise.

Here, we could use ELBO for different number of clones as an indictor for model selection. However, this is still imperfect. One empirical suggestion is to choose the n_clones when ELBO stops increasing dramatically, for example in the case below, we will pick 3 clones.

## Input files
```
CELL_DIR=data/cellSNP_mat
CELL_FILE=data/cells.cellSNP.vcf.gz
DONOR_FILE=data/donors.cellSNP.vcf.gz
DONOR_FILE_PART=data/donors.two.cellSNP.vcf.gz

mkdir data/outs/
```

## MODE 1: no donor genotype
```
OUT_DIR=data/outs/cellSNP_noGT
vireo -c $CELL_DIR -N 4 -o $OUT_DIR --randSeed 2 #--ASEmode # --extraDonor 0
# vireo -c $CELL_DIR -N 4 -o $OUT_DIR --randSeed 2 --extraDonor 1 #--ASEmode
```

## MODE 2: given donor genotype
```
OUT_DIR=data/outs/cellSNP_PL
vireo -c $CELL_FILE -d $DONOR_FILE -o $OUT_DIR -N 4 --randSeed 2 #--genoTag PL
```

## MODE 3: given part donor genotype
```
OUT_DIR=data/outs/cellSNP_part
vireo -c $CELL_FILE -d $DONOR_FILE_PART -o $OUT_DIR -N 4 --randSeed 2
```

## MODE 4: given donor genotype but not perfect, use as prior to learn
```
OUT_DIR=data/outs/cellSNP_learn
vireo -c $CELL_FILE -d $DONOR_FILE -o $OUT_DIR --randSeed 2 -N 4 --forceLearnGT
```

## MODE 5: given donor genotype but too many
```
OUT_DIR=data/outs/cellSNP_PL3
vireo -c $CELL_FILE -d $DONOR_FILE -o $OUT_DIR -N 3 --randSeed 2 #--genoTag PL
```

## Generating genotype barcodes
```
donor_vcf=data/outs/cellSNP_noGT/GT_donors.vireo.vcf.gz
GTbarcode -i $donor_vcf -o data/outs/cellSNP_noGT/GT_barcodes.tsv --randSeed 1 # --noHomoAlt
```