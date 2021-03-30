# Shallow shotgun data pipeline

**>> deprecated, use the metagenomics workflows in single-end mode <<**

**Feasible data:**

- single end metagenomic reads
- low depth (<10M reads per sample)

## Workflow

**Definition**: `main.nf`

### Steps:

### Steps:

1. Adapter and quality trimming with [fastp](https://github.com/OpenGene/fastp)
   *quality reports in HTML and JSON are provided for each file*
2. Co-assembly with [MegaHit](https://github.com/voutcn/megahit)
3. Finding genes *de novo* [prodigal](https://github.com/hyattpd/Prodigal)
4. Transcript-level alignment with [minimap2](https://github.com/lh3/minimap2)
5. Gene quantification (not gene expression!) with the EM counter from [salmon](https://salmon.readthedocs.io/en/latest/salmon.html)
6. Protein annotation using the [EGGNoG mapper](https://github.com/eggnogdb/eggnog-mapper)
7. Replication rates using [iRep](https://www.nature.com/articles/nbt.3704)

