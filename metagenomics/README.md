# Metagenomic workflow

**Feasible data:**

- paired or single end metagenomic shotgun sequencing
- decent depth (>10M reads per sample)

## Basic (functional) workflow

**Definition**: `main.nf`

### Steps:

1. Adapter and quality trimming with [fastp](https://github.com/OpenGene/fastp)
   *quality reports in HTML and JSON are provided for each file*
2. Assembly with [MegaHit](https://github.com/voutcn/megahit)
3. Finding genes *de novo* with [prodigal](https://github.com/hyattpd/Prodigal)
4. Clustering of all genes
5. Pufferfish mapping index creation (needed for next step)
6. Gene quantification (mapping + counting) with [salmon](https://salmon.readthedocs.io/en/latest/salmon.html)
7. Protein annotation using the [EGGNoG mapper](https://github.com/eggnogdb/eggnog-mapper)
8. Generating kmer-based taxonomic profiles using [Kraken2](https://ccb.jhu.edu/software/kraken2/) and [Bracken](https://ccb.jhu.edu/software/bracken/)

Additional workflows can be run *after* the basic workflow has finished.

## Replication rates

**Definition**: replication.nf

1. Alignment to [~3K high quality assemblies](https://www.nature.com/articles/s41586-019-1058-x) from the gut microbiome with bowtie2
2. extraction of coverage maps
3. Quanitifaction of peak-to-trough ratios (PTRs) with [coptr](https://github.com/tyjo/coptr)

## Binning workflow

**Definition**: `binning.nf`


### Steps:

1. metagenomic binning with [Metabat2](https://bitbucket.org/berkeleylab/metabat/)
2. assembly taxonomy assignment using [BAT](https://github.com/dutilh/CAT)
3. quality assessment using [checkM](https://ecogenomics.github.io/CheckM/)

### Setup

```bash
conda env create -f conda.yml
```

## Options:

```
~~~ Gibbons Lab Metagenomics Workflow ~~~

Usage:
A run using all,default parameters can be started with:
> nextflow run main.nf --resume

A run with all parametrs set would look like:
> nextflow run main.nf --data_dir=./data --refs=/my/references --single_end=false \
                       --trim_front=5 --min_length=50 --quality_threshold=20 --read_length=150 --threshold=10

General options:
  --data_dir [str]              The main data directory for the analysis (must contain `raw`).
  --read_length [str]           The length of the reads.
  --single_end [bool]           Specifies that the input is single-end reads.
Reference DBs:
  --refs [str]                  Folder in which to find references DBs.
  --eggnogg_refs [str]          Where to find EGGNOG references. Defaults to <refs>/eggnog.
  --kraken2_db [str]            Where to find the Kraken2 reference. Defaults to <refs>/kraken2_default.
Quality filter:
  --trim_front [str]            How many bases to trim from the 5' end of each read.
  --min_length [str]            Minimum accepted length for a read.
  --quality_threshold [str]     Smallest acceptable average quality.
  --threshold [str]             Smallest abundance threshold (in reads) used by Kraken.
```

