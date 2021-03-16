#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${baseDir}/data"
params.single_end = false
params.min_contig_length = 5000

maxcpus = 24

process bin_align {
    cpus maxcpus

    input:
    tuple val(id), file(reads)
    path(assembly)

    output:
    path("${id}.bam")

    """
    minimap2 -ax sr -N 100 -t ${task.cpus} ${assembly} ${reads} | \
    samtools sort -@${task.cpus} -o ${id}.bam && \
    samtools index ${id}.bam ${id}.bai
    """
}

process bin {
    cpus maxcpus
    publishDir "${params.data_dir}"

    input:
    path(bam_files)
    path(assembly)

    output:
    path("metabat2_bins")

    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bam_files} && \
    metabat2 -i ${assembly} -a depth.txt -o metabat2_bins/bin \
        -t ${task.cpus} -m 2500 -s 100000
    """
}

process build_database {
    cpus maxcpus
    publishDir "${params.data_dir}/refs"

    output:
    tuple path("CAT_database"), path("CAT_taxonomy")

    """
    CAT prepare --fresh -d CAT_database -t CAT_taxonomy -n ${task.cpus}
    """
}

process filter_contigs {
    cpus 1
    publishDir "${params.data_dir}/assembled/contigs"

    input:
    path(contigs)

    output:
    path("long_contigs.fna")

    """
    #!/usr/bin/env Rscript

    library(Biostrings)

    contigs <- readDNAStringSet("${contigs}")
    writeXStringSet(contigs[width(contigs) > ${params.min_contig_length}], "long_contigs.fna")
    """
}

process bin_taxonomy {
    cpus maxcpus
    publishDir "${params.data_dir}"

    input:
    path(bins)
    tuple path(database_folder), path(taxonomy_folder)

    output:
    path("bins.classification.names.txt")

    """
    CAT bins -r 10 -f 0.1 -b ${bins} -s .fa -d ${database_folder} \
        -t ${taxonomy_folder} -o bins -n ${task.cpus} && \
    CAT add_names -i bins.bin2classification.txt -o bins.classification.names.txt -t ${taxonomy_folder}
    """
}

process contig_taxonomy {
    cpus maxcpus
    publishDir "${params.data_dir}"

    input:
    path(contigs)
    tuple path(database_folder), path(taxonomy_folder)

    output:
    tuple path("contigs.classification.names.txt"), path("contigs.summary.txt")

    """
    CAT contigs -c ${contigs} -d ${database_folder} -t ${taxonomy_folder} \
        -n ${task.cpus} --out_prefix contigs && \
    CAT add_names -i contigs.contig2classification.txt -o contigs.classification.names.txt \
        -t ${taxonomy_folder} --only_official && \
    CAT summarise -c ${contigs} -i contigs.classification.names.txt -o contigs.summary.txt
    """
}

workflow {
    if (params.single_end) {
        Channel
            .fromPath("${params.data_dir}/preprocessed/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{reads}
    } else {
        Channel
            .fromFilePairs([
                "${params.data_dir}/preprocessed/*_R{1,2}_001.fastq.gz",
                "${params.data_dir}/preprocessed/*_{1,2}.fastq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${params.data_dir}!" }
            .set{reads}
    }

    Channel.from("${params.data_dir}/assembled/contigs/final.contigs.fa").set{assembly}

    dbs = build_database()

    aligned = bin_align(reads, assembly)
    binned = bin(aligned.collect(), assembly)
    bin_classified = bin_taxonomy(binned, dbs)
    contig_classified = contig_taxonomy(filter_contigs(assembly), dbs)
}
