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

For CARVEME use conda:

```bash
conda env create -f conda.yml
```

For gapseq use singularity:

```bash
nextflow run main.nf --method gapseq -with-singularity cdiener/gapseq:2021.04.5
```

Make sure that singularity runs with the `--no-home` option and map your data dir
into singularity. For instance, by adding the following to your nextflow config:

```
singularity.runOptions = "--no-home -B /MY/DATA/DIR,/tmp"
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
  --method                      The method used to build the models. Either `carveme` or `gapseq`.
Growth Media:
  --media_db                    A file containing growth media specification for CARVEME. This is a
                                TSV file for carveme and a CSV file for gapseq.
  --media                       Comma-separated list of media names to use. Only used for CARVEME.
```
