# Metabolic Model Builder

**Feasible data:**

- good quality genome assemblies in FASTA format

## Basic (functional) workflow

**Definition**: `main.nf`

### Steps:

1. Initialize the CARVEME DB if necessary
2. Finding genes *de novo* [prodigal](https://github.com/hyattpd/Prodigal)
3. Model construction using [CARVEME](https://carveme.readthedocs.io/)
4. Model quality checks using [MEMOTE](https://memote.readthedocs.io/))


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
