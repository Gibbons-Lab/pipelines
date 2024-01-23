#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data = "data"
params.read_length = 150
params.single_end = false
params.min_contig_len = 1000
params.min_length = 60
params.min_quality = 20
params.trim_front = 5
params.threads = 24
params.version = 1
params.gtdb = "/proj/gibbons/refs/gtdbtk_r207v2"
params.checkm2 = "/proj/gibbons/refs/checkm2/uniref100.KO.1.dmnd"

def helpMessage() {
    log.info"""
    ~~~ Gibbons Lab Isolate Shotgun Sequencing Workflow ~~~

    Usage:
    A run can be started with:
    > nextflow run main.nf.nf --resume --reference e_coli_k12.fna

    General options:
      --data [str]                  The main data directory for the analysis (must contain `raw`).
      --read_length [str]           The length of the reads.
      --single_end [bool]           Specifies that the input is single-end reads.
      --threads [int]               The maximum number of threads a single process can use.
                                    This is not the same as the maximum number of total threads used.
    Quality filter:
      --trim_front [str]            How many bases to trim from the 5' end of each read.
      --min_length [str]            Minimum accepted length for a read.
      --quality_threshold [str]     Smallest acceptable average quality.
    Assembly:
      --min_contig_len [int]        Minimum length of a contig.
      --mincov                      Minimum length coverage required.
      --minsamples [int]            Minimum number of samples with good mapping.
    Taxonomy:
      --gtdb                        Path to the downloaded GTDB-TK files.
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
    memory "4 GB"
    publishDir "${baseDir}/${params.data}/preprocessed"

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
            -3 -M ${params.min_quality} -r -w ${task.cpus} \
            --max_len1 ${params.read_length}
        """

    else
        """
        fastp -i ${reads[0]} -I ${reads[1]} \
            -o ${id}_filtered_R1.fastq.gz -O ${id}_filtered_R2.fastq.gz\
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.min_quality} -r -w ${task.cpus} \
            --max_len1 ${params.read_length} --max_len2 ${params.read_length}
        """
}

process megahit {
    cpus 4
    memory "32 GB"
    publishDir "${baseDir}/${params.data}/assemblies", mode: "copy", overwrite: true

    input:
    tuple val(id), path(reads), path(json), path(report)

    output:
    tuple val(id), path("${id}.asm1.fna")

    script:
    if (params.single_end)
        """
        megahit -r ${reads} -o contigs -t ${task.cpus} -m 0.5 \
                --min-contig-len ${params.min_contig_len} --out-prefix ${id}
        sed -i -e "s/^>/>${id}_/" contigs/${id}.contigs.fa
        mv contigs/${id}.contigs.fa ${id}.asm${params.version}.fna
        """
    else
        """
        megahit -1 ${reads[0]} -2 ${reads[1]} -o contigs -t ${task.cpus} -m 0.5 \
                --min-contig-len ${params.min_contig_len} --out-prefix ${id}
        sed -i -e "s/^>/>${id}_/" contigs/${id}.contigs.fa
        mv contigs/${id}.contigs.fa ${id}.asm${params.version}.fna
        """
}

process multiqc {
    cpus 1
    memory "8 GB"
    publishDir "${baseDir}/${params.data}", mode: "copy", overwrite: true

    input:
    val(preprocess)

    output:
    path("multiqc_report.html")

    """
    multiqc ${baseDir}/data/preprocessed
    """
}

process classify {
    cpus params.threads
    memory "64 GB"
    publishDir "${baseDir}/${params.data}", mode: "copy", overwrite: true

    input:
    path(assemblies)

    output:
    path("gtdb")

    """
    mkdir assemblies && mv *.asm1.fna assemblies
    GTDBTK_DATA_PATH=${params.gtdb} gtdbtk classify_wf \
        --genome_dir assemblies/ --prefix assemblies \
        --cpus ${task.cpus} --out_dir gtdb
    """
}

process tree {
    cpus params.threads
    memory "64 GB"
    publishDir "${baseDir}/${params.data}", mode: "copy", overwrite: true

    input:
    path(gtdb)

    output:
    path("isolates.nwk")

    """
    gunzip -k ${gtdb}/align/assemblies.bac120.user_msa.fasta.gz && \
        mv ${gtdb}/align/assemblies.bac120.user_msa.fasta msa.fna
    OMP_NUM_THREADS=${task.cpus} FastTreeMP \
        -wag -out isolates.nwk msa.fna && \
        rm msa.fna
    """
}

process checkm {
    cpus params.threads
    memory "64 GB"
    publishDir "${baseDir}/${params.data}", mode: "copy", overwrite: true

    input:
    path(assemblies)

    output:
    path("checkm")

    """
    checkm2 predict \
        --database_path ${params.checkm2} \
        --threads ${task.cpus} \
        --input ${assemblies} \
        --output-directory checkm
    """
}

workflow {
    data_dir = "${baseDir}/${params.data}"
    // find files
    if (params.single_end) {
        Channel
            .fromPath("${data_dir}/raw/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{raw}
    } else {
        Channel
            .fromFilePairs([
                "${data_dir}/raw/*_R{1,2}_001.fastq.gz",
                "${data_dir}/raw/*_{1,2}.fastq.gz",
                "${data_dir}/raw/*_R{1,2}.fastq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${data_dir}/raw!" }
            .set{raw}
    }

    raw | preprocess | megahit
    multiqc(preprocess.out.collect())

    megahit.out
        .map{it -> it[1]}
        .filter{it.size() > 0}
        .collect()
        .set{assemblies}
    assemblies | checkm
    assemblies| classify | tree
}
