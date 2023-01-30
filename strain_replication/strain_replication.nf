nextflow.enable.dsl=2

params.reference = "reference.fna"
params.decoy = "decoy.fna"
params.single_end = false
params.trim_front = 5
params.min_length = 50
params.quality_threshold = 20
params.threads = 12
params.minreads = 5000
params.mincov = 0.75
params.minsamples = 1

def helpMessage() {
    log.info"""
    ~~~ Gibbons Lab Strain Replication Workflow ~~~

    Usage:
    A run using all,default parameters can be started with:
    > nextflow run strain_replication.nf --resume --reference e_coli_k12.fna

    General options:
      --data_dir [str]              The main data directory for the analysis (must contain `raw`).
      --read_length [str]           The length of the reads.
      --single_end [bool]           Specifies that the input is single-end reads.
      --threads [int]               The maximum number of threads a single process can use.
                                    This is not the same as the maximum number of total threads used.
    Quality filter:
      --trim_front [str]            How many bases to trim from the 5' end of each read.
      --min_length [str]            Minimum accepted length for a read.
      --quality_threshold [str]     Smallest acceptable average quality.
      --threshold [str]             Smallest abundance threshold used by Kraken.
    COPTR:
      --minreads [int]              Minimum number of reads covering a genome.
      --mincov                      Minimum length coverage required.
      --minsamples [int]            Minimum number of samples with good mapping.
    """.stripIndent()
}

params.help = false
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

process preprocess {
    cpus 4
    publishDir "${baseDir}/data/preprocessed"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_filtered_R*.fastq.gz"), path("${id}_fastp.json"), path("${id}.html")

    script:
    if (params.single_end)
        """
        fastp -i ${reads[0]} -o ${id}_filtered_R1.fastq.gz \
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r -w ${task.cpus}
        """

    else
        """
        fastp -i ${reads[0]} -I ${reads[1]} \
            -o ${id}_filtered_R1.fastq.gz -O ${id}_filtered_R2.fastq.gz\
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r -w ${task.cpus}
        """
}

process build_index {
    cpus params.threads
    publishDir "${baseDir}/data/"

    input:
    path(ref)
    path(decoy)

    output:
    path("reference_index")

    """
    mkdir genomes && cp ${ref} ${decoy} genomes
    mkdir reference_index
    coptr index --bt2-threads ${task.cpus} genomes reference_index/index
    """
}

process map_reads {
    cpus params.threads/4
    publishDir "${baseDir}/data/mapped"

    input:
    tuple val(id), path(reads), path(report), path(html)
    each path(index)

    output:
    tuple val(id), path("${id}.bam")

    script:
    if (params.single_end)
        """
      mkdir files bam
      (cd files && ln -s ../${reads[0]} reads_1.fastq.gz)
      coptr map --threads ${task.cpus} ${index}/index files bam
      mv bam/*.bam ${id}.bam
      """
    else
      """
      mkdir files bam
      (cd files && ln -s ../${reads[0]} reads_1.fastq.gz && ln -s ../${reads[1]} reads_2.fastq.gz)
      coptr map --threads ${task.cpus} --paired ${index}/index files bam
      mv bam/*.bam ${id}.bam
      """
}

process multiqc {
    publishDir "${baseDir}/data", mode: "copy", overwrite: true

    input:
    path(files)

    output:
    path("multiqc_report.html")

    """
    multiqc .
    """
}

process alignment_stats {
  cpus 1
  publishDir "${baseDir}/data/stats"

  input:
  tuple val(id), path(bam)

  output:
  path("${id}.sts")

  """
  samtools stats ${bam} > ${id}.sts
  """

}

process extract_coverage {
  cpus 1
  publishDir "${baseDir}/data"

  input:
  tuple val(id), path(bam)

  output:
  path("coverage/*.cm.pkl")

  """
  mkdir coverage
  coptr extract . coverage
  """
}

process estimate_ptr {
  cpus 1
  publishDir "${baseDir}/data",  mode: "copy", overwrite: true

  input:
  path(coverage)

  output:
  path("rates.csv")

  """
  coptr estimate \\
    --min-reads ${params.minreads} \\
    --min-cov ${params.mincov} \\
    --min-samples ${params.minsamples} . rates.csv
  """
}

process count_reads {
  cpus 1
  publishDir "${baseDir}/data",  mode: "copy", overwrite: true

  input:
  path(coverage)

  output:
  path("reference_counts.csv")

  """
  #!/usr/bin/env python

  import pickle
  import pandas as pd

  files = "${coverage}".split()
  counts = []
  for f in files:
    print(f"Processing coverage profile {f}...")
    cov = pickle.load(open(f, "rb"))
    c = pd.Series({id: c.count_reads() for id, c in cov.items()})
    c["sample_id"] = f.split(".cm.pkl")[0]
    counts.append(c)
  counts = pd.DataFrame.from_records(counts)
  counts.to_csv("reference_counts.csv", index=False)
  """
}

workflow {
// find files
    if (params.single_end) {
        Channel
            .fromPath("${baseDir}/data/raw/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{reads}
    } else {
        Channel
            .fromFilePairs([
                "${baseDir}/data/raw/*_filtered_R{1,2}.fastq.gz",
                "${baseDir}/data/raw/*_R{1,2}_001.fastq.gz",
                "${baseDir}/data/raw/*_{1,2}.fastq.gz",
                "${baseDir}/data/raw/*_{1,2}.fq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${baseDir}/data!" }
            .set{reads}
    }

    Channel
        .fromPath("${baseDir}/data/${params.reference}")
        .ifEmpty { error "Cannot find the reference in ${baseDir}/data!" }
        .set{ref}
    Channel
        .fromPath("${baseDir}/data/${params.decoy}")
        .ifEmpty { error "Cannot find the decoy in ${baseDir}/data!" }
        .set{decoy}
    build_index(ref, decoy)
    reads | preprocess
    map_reads(preprocess.out, build_index.out) | alignment_stats

    fastp_reps = preprocess.out.map{it -> it[2]}.collect()
    alignment_reps = alignment_stats.out.collect()

    multiqc(fastp_reps.concat(alignment_reps).collect())

    map_reads.out | extract_coverage
    estimate_ptr(extract_coverage.out.collect())
    extract_coverage.out.collect() | count_reads
}
