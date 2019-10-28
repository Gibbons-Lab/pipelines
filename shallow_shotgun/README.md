## Shallow shotgun data pipeline

**Feasible data:**

- single end metagenomic reads
- low depth (<10M reads per sample)

### Steps:

1. Adapter and quality trimming with [fastp](https://github.com/OpenGene/fastp)
   *quality reports in HTML and JSON are provided for each file*
2. Abundance quantification with [Kraken2](https://ccb.jhu.edu/software/kraken2/) + [Bracken](https://ccb.jhu.edu/software/bracken/)
3. Merging of rank-level abundance tables across samples
