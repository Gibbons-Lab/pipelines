#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${baseDir}/data"
params.raw_data = "raw"
params.transcripts = "data/transcripts.fna.gz"

params.single_end = false
params.trim_front_fwd = 5
params.trim_front_rev = 5
params.min_length = 50
params.quality_threshold = 20
params.read_length = 150
params.threads = 20

def helpMessage() {
    log.info"""
    ~~~ Gibbons Lab Metagenomics Workflow ~~~

    Usage:
    A run using all,default parameters can be started with:
    > nextflow run main.nf --resume

    A run with all parametrs set would look like:
    > nextflow run main.nf --data_dir=./data --single_end=false --refs=/my/references --single_end=false \\
                           --trim_front=5 --min_length=50 --quality_threshold=20 --read_length=150 --threshold=10

    General options:
      --data_dir [str]              The main data directory for the analysis (must contain `raw`).
      --read_length [str]           The length of the reads.
      --single_end [bool]           Specifies that the input is single-end reads.
      --threads [int]               The maximum number of threads a single process can use.
                                    This is not the same as the maximum number of total threads used.
    Reference DBs:
      --refs [str]                  Folder in which to find references DBs.
      --eggnogg_refs [str]          Where to find EGGNOG references. Defaults to <refs>/eggnog.
      --kraken2_db [str]            Where to find the Kraken2 reference. Defaults to <refs>/kraken2_default.
    Quality filter:
      --trim_front [str]            How many bases to trim from the 5' end of each read.
      --min_length [str]            Minimum accepted length for a read.
      --quality_threshold [str]     Smallest acceptable average quality.
      --threshold [str]             Smallest abundance threshold used by Kraken.
    Assembly:
      --contig_length [int]         Minimum length of a contig.
      --identity [double]           Minimum average nucleotide identity.
      --overlap [double]            Minimum required overlap between contigs.
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
    publishDir "${params.data_dir}/preprocessed"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_filtered_R*.fastq.gz"), path("${id}_fastp.json"), path("${id}.html")

    script:
    if (params.single_end)
        """
        fastp -i ${reads[0]} -o ${id}_filtered_R1.fastq.gz \
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front_fwd} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r -w ${task.cpus} \
            --max_len1 ${params.read_length}
        """

    else
        """
        fastp -i ${reads[0]} -I ${reads[1]} \
            -o ${id}_filtered_R1.fastq.gz -O ${id}_filtered_R2.fastq.gz\
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front_fwd} --trim_front2 ${params.trim_front_rev} \
            -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r -w ${task.cpus} \
            --max_len1 ${params.read_length} --max_len2 ${params.read_length}
        """
}

process multiqc {
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(preprocessed)
    path(salmon)

    output:
    path("multiqc_report.html")

    """
    multiqc ${params.data_dir}/preprocessed ${params.data_dir}/salmon
    """
}


process index {
    cpus params.threads
    publishDir "${projectDir}/data"

    input:
    path(transcripts)

    output:
    path("salmon_index")

    """
    salmon index -p ${task.cpus} -t ${transcripts} -i salmon_index
    """
}

process quantify {
    cpus 4
    memory "64 GB"

    publishDir "${projectDir}/data/salmon"

    input:
    tuple val(id), path(reads), path(json), path(html)
    each path(index)

    output:
    path("${id}")

    script:
    if (params.single_end)
        """
        salmon quant --meta -p ${task.cpus} -l A -i ${index} -r ${reads} -o ${id} --validateMappings
        """
    else
        """
        salmon quant --meta -p ${task.cpus} -l A -i ${index} -1 ${reads[0]} -2 ${reads[1]} -o ${id} --validateMappings
        """
}

process merge_counts {
    cpus 1
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(salmon_quants)

    output:
    path("transcript_counts.csv.gz")

    """
    #!/usr/bin/env python

    from sys import stdin
    from os import path
    import pandas as pd
    import gzip

    paths = "${salmon_quants}"
    paths = [path.join(p, "quant.sf") for p in paths.split(" ")]
    nread = 0
    with gzip.open("transcript_counts.csv.gz", "ab") as gzf:
        for p in paths:
            sample = path.splitext(path.basename(p))[0]
            print("Processing sample {sample}...")
            try:
                counts = pd.read_csv(p, sep="\t").query("NumReads > 0.1")
            except Exception:
                continue
            nread += 1

            counts.columns = [
                "locus_tag", "length", "effective_length", "tpm", "reads"]
            counts["sample_id"] = sample
            print(f"writing compressed output for sample {sample}...")
            counts.to_csv(gzf, header=(nread==1),
                          index=False)
    """
}

workflow {
    // find files
    if (params.single_end) {
        Channel
            .fromPath("${params.data_dir}/${params.raw_data}/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{raw}
    } else {
        Channel
            .fromFilePairs([
                "${params.data_dir}/raw/*_R{1,2}_001.fastq.gz",
                "${params.data_dir}/raw/*_{1,2}.fastq.gz",
                "${params.data_dir}/raw/*_{1,2}.fq.gz",
                "${params.data_dir}/raw/*_R{1,2}.fastq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${params.data_dir}/${params.raw_data}!" }
            .set{raw}
    }

    // quality filtering
    preprocess(raw)

    // build the Salmon index
    Channel.fromPath("${params.transcripts}") | index

    // Quantify the transcripts
    quantify(preprocess.out, index.out) | merge_counts

    // quality overview
    multiqc(preprocess.out.map{it[2]}.collect(), quantify.out.collect())
}