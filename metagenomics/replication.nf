#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${baseDir}/data"
params.single_end = false
params.min_reads = 5000
params.IGG = "/proj/gibbons/refs/IGG_v1.0_split"
params.threads = 12

def helpMessage() {
    log.info"""
    ~~~ Gibbons Lab Replication Rate workflow ~~~

    Usage:
    A run using all,default parameters can be started with:
    > nextflow run replication.nf --resume

    A run with all parametrs set would look like:
    > nextflow run main.nf --data_dir=./data --single_end true --threads 12 --min_reads 2500 --IGG /refs/IGG

    General options:
      --data_dir [str]              The main data directory for the analysis (must contain `raw`).
      --single_end [bool]           Whether the data is single-end sequencing data.
      --threads [int]               The maximum number of threads a single process can use.
                                    This is not the same as the maximum number of total threads used.

    COPTR options:
      --min_reads [int]             Minimum number of reads for a genome to calculate PTRs.
      --IGG [path]                  Location of the IGG bowtie2 reference.
    """.stripIndent()
}

params.help = false
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

process map_reads {
    cpus 4
    publishDir "${params.data_dir}/alignments"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.bam")

    script:
    if (params.single_end)
        """
      mkdir files bam1 bam2
      (cd files && ln -s ../${reads[0]} reads_1.fastq.gz)
      coptr map --threads ${task.cpus} ${params.IGG}/IGG_v1.0-1 files bam1
      coptr map --threads ${task.cpus} ${params.IGG}/IGG_v1.0-2 files bam2
      coptr merge bam1/*.bam bam2/*.bam ${id}.bam
      """
    else
      """
      mkdir files bam1 bam2
      (cd files && ln -s ../${reads[0]} reads_1.fastq.gz && ln -s ../${reads[1]} reads_2.fastq.gz)
      coptr map --threads ${task.cpus} --paired ${params.IGG}/IGG_v1.0-1 files bam1
      coptr map --threads ${task.cpus} --paired ${params.IGG}/IGG_v1.0-2 files bam2
      coptr merge bam1/*.bam bam2/*.bam ${id}.bam
      """
}

process extract_coverage {
  cpus 1
  publishDir "${params.data_dir}/coverage"

  input:
  tuple val(id), path(bam)

  output:
  path("coverage/*.*")

  """
  mkdir coverage
  coptr extract . coverage
  """
}

process estimate_ptr {
  cpus 1
  publishDir "${params.data_dir}",  mode: "copy", overwrite: true

  input:
  path(coverage)

  output:
  path("rates.csv")

  """
  coptr estimate --min-reads ${params.min_reads} . rates.csv
  """
}

process annotate_ptr {
  cpus 1
  publishDir "${params.data_dir}",  mode: "copy", overwrite: true

  input:
  path(rates)

  output:
  path("annotated_rates.csv")

  shell:
  '''
  #!/usr/bin/env python

  import pandas as pd

  rates = pd.read_csv("!{rates}")
  meta = pd.read_csv("https://portal.nersc.gov/cfs/m342/HGM/IGG_species_info_23790.tsv", sep="\t")
  rates.rename(columns={"log2(PTR):genome_id/sample_id": "genome_id"}, inplace=True)
  rates = rates.melt(id_vars="genome_id", var_name="sample_id", value_name="log2_ptr").dropna()
  rates["representative_genome"] = rates.genome_id.str.replace("\\.\\w+$", "", regex=True)

  taxa = meta.gtdb_taxonomy.str.replace("\\w__", "", regex=True).str.split(";", expand=True)
  taxa.columns = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
  meta = pd.concat([meta, taxa], axis=1)

  merged = pd.merge(rates, meta, on="representative_genome")
  merged.to_csv("annotated_rates.csv", index=False)
  '''
}

workflow {
  // find files
    if (params.single_end) {
        Channel
            .fromPath("${params.data_dir}/preprocessed/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{reads}
    } else {
        Channel
            .fromFilePairs([
                "${params.data_dir}/preprocessed/*_filtered_R{1,2}.fastq.gz",
                "${params.data_dir}/preprocessed/*_R{1,2}_001.fastq.gz",
                "${params.data_dir}/preprocessed/*_{1,2}.fastq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${params.data_dir}!" }
            .set{reads}
    }

    reads | map_reads | extract_coverage
    estimate_ptr(extract_coverage.out.collect()) | annotate_ptr
}
