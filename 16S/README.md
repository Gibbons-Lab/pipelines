# 16S amplicon sequencing workflow

**Feasible data:**

- single-end or paired-end 16S amplicon sequencing data
- decent depth (>10K reads per sample)

## Basic workflow

**Definition**: `main.nf`

### Steps:

1. Automatic detection of Illumina read files and management of multiple runs
2. Trimming, filtering and quality metrics (base qualities, entropy, lengths) with DADA2 and mbtools
3. Denoising using DADA2
4. Taxonomy assignemnt using DADA2 and SILVA
5. Creation of tabular and phyloseq output
