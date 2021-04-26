#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${baseDir}/data"
params.single_end = false
params.min_reads = 5000
params.IGG = "/proj/gibbons/refs/IGG_v1.0_split"

def max_threads = Runtime.runtime.availableProcessors()

def helpMessage() {
    log.info"""
    ~~~ Gibbons Lab Replication Rate workflow ~~~

    Usage:
    A run using all,default parameters can be started with:
    > nextflow run replication.nf --resume

    A run with all parametrs set would look like:
    > nextflow run main.nf --data_dir=./data --media_db=media.tsv --media="LB,M9"

    General options:
      --data_dir [str]              The main data directory for the analysis (must contain `raw`).
      --single_end [bool]           Whether the data is single-end sequencing data.

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

process download_coptr {
    cpus 1

    output:
    path("coptr")

    """
    wget https://github.com/tyjo/coptr/archive/refs/heads/master.zip
    unzip master.zip
    mv coptr-master coptr
    """
}

process map_reads {
    cpus 4
    publishDir "${params.data_dir}/alignments"

    input:
    tuple val(id), path(reads), path(coptr)

    output:
    tuple val(id), path("${id}.bam")

    script:
    if (params.single_end)
        """
      mkdir files bam1 bam2
      (cd files && ln -s ../${reads[0]} ${id}_1.fastq.gz)
      python ${coptr}/coptr.py map --threads ${task.cpus} ${params.IGG}/IGG_v1.0-1 files bam1
      python ${coptr}/coptr.py map --threads ${task.cpus} ${params.IGG}/IGG_v1.0-2 files bam2
      python ${coptr}/coptr.py merge bam1/*.bam bam2/*.bam ${id}.bam
      """
    else
      """
      mkdir files bam1 bam2
      (cd files && ln -s ../${reads[0]} ${id}_1.fastq.gz && ln -s ../${reads[1]} ${id}_2.fastq.gz)
      python ${coptr}/coptr.py map --threads ${task.cpus} --paired ${params.IGG}/IGG_v1.0-1 files bam1
      python ${coptr}/coptr.py map --threads ${task.cpus} --paired ${params.IGG}/IGG_v1.0-2 files bam2
      python ${coptr}/coptr.py merge bam1/*.bam bam2/*.bam ${id}.bam
      """
}

process extract_coverage {
  cpus 1
  publishDir "${params.data_dir}/coverage"

  input:
  tuple val(id), path(bam), path(coptr)

  output:
  path("coverage/${id}.cm.pkl")

  """
  mkdir coverage
  python ${coptr}/coptr.py extract . coveraged
  """
}

process estimate_ptr {
  cpus max_threads
  publishDir "${params.data_dir}",  mode: "copy", overwrite: true

  input:
  tuple path(coverage), path(coptr)

  output:
  path("rates.csv")

  """
  python ${coptr}/coptr.py estimate . rates.csv --min-reads ${params.min_reads}
  """
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

    download_coptr()
    reads.combine(download_coptr.out) | map_reads
    map_reads.out.combine(download_coptr.out) | extract_coverage
    extract_coverage.out.collect().combine(download_coptr.out) | estimate_ptr
}
