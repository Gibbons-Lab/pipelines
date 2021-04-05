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
~~~ Gibbons Lab Metabolic Model Builder Workflow ~~~

Usage:
A run using all,default parameters can be started with:
> nextflow run main.nf --resume

A run with all parametrs set would look like:
> nextflow run main.nf --data_dir=./data --media_db=media.tsv --media="LB,M9"

General options:
  --data_dir [str]              The main data directory for the analysis (must contain `raw`).
Growth Media:
  --media_db                    A file containing growth media specification for CARVEME.
  --media                       Comma-separated list of media names to use.
```
