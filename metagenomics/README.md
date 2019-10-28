## Metagenomic workflow

**Feasible data:**

- paired end metagenomic shotgun sequencing
- decent depth (>10M reads per sample)

### Steps:

1. Adapter and quality trimming with [fastp](https://github.com/OpenGene/fastp)
   *quality reports in HTML and JSON are provided for each file*
2. Co-assembly with [MegaHit](https://github.com/voutcn/megahit)
3. Annotation with [PROKKA](https://github.com/tseemann/prokka)
4. Transcript-level alignment with [minimap2](https://github.com/lh3/minimap2)
5. Gene quantification (not gene expression!) with the EM counter from [salmon](https://salmon.readthedocs.io/en/latest/salmon.html)

Optional steps:

- metagenomic binning with [Metabat2](https://bitbucket.org/berkeleylab/metabat/)
- binning quality check with [checkm](https://ecogenomics.github.io/CheckM/)
