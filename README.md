# :hammer: :triangular_ruler: Pipelines

This repo contains various analysis pipelines for the lab. Here are the basic
rules:

- each folder includes pipelines for a particular analysis - data type combination
- pipelines are [nextflow](https://www.nextflow.io/) workflows
- each pipeline comes with a list of conda environment files that manage the required software

## Data layout

Pipelines will usually operate from a top level project
directory structured in the following way:

```
project/
    > [WORKFLOW].nf
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
conda env create -f eggnog.yml
```

Activate the environment (usualy named after the pipeline):

```bash
conda activate metagenomics
```

If present also install R dependencies

```bash
Rscript setup.R
```

You may also create a [nextflow config](https://www.nextflow.io/docs/latest/config.html) either in the project
directory as `nextflow.config` or in your user HOME as `~/.nextflow/config`. Here is an example for moneta

```
executor {
    name = "local"
    cpus = 32
    memory = "350GB"
}

tower {
    accessToken = "<TOKEN>"
    enabled = false
}
```

The tower part is only necessary if you want to use [Nextflow Tower](https://tower.nf)
to track your pipeline.

## Run the pipeline

After setup you can run the pipeline with

```bash
nextflow run [WORKFLOW].nf -resume
```

By default this will use all available CPUs and RAM unless specified otherwise in a personal [netxflow config](https://www.nextflow.io/docs/latest/config.html#scope-executor).
