# Metagenomic workflow

**Feasible data:**

- paired or single end metagenomic shotgun sequencing
- decent depth (>10M reads per sample)

## Basic (functional) workflow

**Definition**: `main.nf`

### Steps:

1. Adapter and quality trimming with [fastp](https://github.com/OpenGene/fastp)
   *quality reports in HTML and JSON are provided for each file*
2. Co-assembly with [MegaHit](https://github.com/voutcn/megahit)
3. Finding genes *de novo* [prodigal](https://github.com/hyattpd/Prodigal)
4. Transcript-level alignment with [minimap2](https://github.com/lh3/minimap2)
5. Gene quantification (not gene expression!) with the EM counter from [salmon](https://salmon.readthedocs.io/en/latest/salmon.html)
6. Protein annotation using the [EGGNoG mapper](https://github.com/eggnogdb/eggnog-mapper)
7. Replication rates using [iRep](https://www.nature.com/articles/nbt.3704)

## Binning workflow

**Definition**: `binning.nf`


### Steps:

- metagenomic binning with [Metabat2](https://bitbucket.org/berkeleylab/metabat/)
- contig taxonomy assignment using [CAT](https://github.com/dutilh/CAT)
- binning quality check with [BAT](https://github.com/dutilh/CAT)

### Setup

```bash
conda env create -f conda.yml
conda env create -f eggnog

Rscript setup.R
```

## Options:

```
~~~ Gibbons Lab Metagenomics Workflow ~~~

Usage:
A run using all,default parameters can be started with:
> nextflow run main.nf --resume

A run with all parametrs set would look like:
> nextflow run main.nf --data_dir=./data --single_end=false --refs=/my/references --single_end=false \
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
  --threshold [str]             Smallest abundance threshold used by Kraken.
```
