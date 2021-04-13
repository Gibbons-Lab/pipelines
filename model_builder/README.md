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
