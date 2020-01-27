# Shallow shotgun data pipeline

**Feasible data:**

- single end metagenomic reads
- low depth (<10M reads per sample)

## Abundance workflow

**Definition**: `Snakefile`

### Steps:

1. Adapter and quality trimming with [fastp](https://github.com/OpenGene/fastp)
   *quality reports in HTML and JSON are provided for each file*
2. Abundance quantification with [Kraken2](https://ccb.jhu.edu/software/kraken2/) + [Bracken](https://ccb.jhu.edu/software/bracken/)
3. Merging of rank-level abundance tables across samples

## De novo workflow

**Definition**: `denovo.nf`

### Steps:

1. Adapter and quality trimming with [fastp](https://github.com/OpenGene/fastp)
   *quality reports in HTML and JSON are provided for each file*
2. Co-assembly with [MegaHit](https://github.com/voutcn/megahit)
3. Finding genes *de novo* [prodigal](https://github.com/hyattpd/Prodigal)
4. Transcript-level alignment with [minimap2](https://github.com/lh3/minimap2)
5. Gene quantification (not gene expression!) with the EM counter from [salmon](https://salmon.readthedocs.io/en/latest/salmon.html)
6. Protein annotation using the [EGGNoG mapper](https://github.com/eggnogdb/eggnog-mapper)
