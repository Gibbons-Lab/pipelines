# Single strain replication rates

**Feasible data:**

- paired or single end isolate shotgun sequencing data at various growth phases
- decent depth (>1M reads per sample)
- provided a reference genome or assembly

## Basic workflow

**Definition**: `strain_replication.nf`

### Steps:

1. Adapter and quality trimming with [fastp](https://github.com/OpenGene/fastp)
   *quality reports in HTML and JSON are provided for each file*
2. Bowtie index build for reference genome
3. Mapping of reads to index
4. Extraction of coverage profiles
5. Estimation of PTRs with [COPTR](https://github.com/tyjo/coptr)
