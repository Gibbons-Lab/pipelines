# :hammer: :triangular_ruler: Pipelines

This repo contains various analysis pipelines for the lab. Here are the basic
rules:

- each folder is a pipeline for a particular analysis - data type combination
- pipelines are [snakemake](https://snakemake.readthedocs.io/en/stable/) workflows
- each pipeline comes with a list of conda environment files that manage the required software
- each pipeline comes with a config file setting some global options

## Data layout

Pipelines will usually operate from a top level project
directory structured in the following way:

`
project/
    > Snakefile
    > config.yml
    > data
        > raw
        > step1
        > step2
        > output1.csv
        > ...
    > figures
        > figure1.png
        > ...
```

The initial raw data lives in `data/raw` and all analysis artifacts should
be written into `data/` as well. Figures go into `figures/`.

## Setup

The first step is to copy or symlink the pipeline files into the top project
directory. After that you can set up a conda environment that includes all software
for the pipeline.

```bash
conda env create -f conda.yml
```

Some pipelines may have separated conda environments for individual steps.
Those are managed by snakemake but will require to run snakemake with the
`--use-conda` flag.

## Run the pipeline

After adjusting values in `config.yml` you can run the pipeline with

```bash
snakemake
```

In general you will want to specify the number of threads used by snakemake. Snakemake
will then balance individual steps over all threads automatically. For instance:

```bash
snakemake --cores 24
```
