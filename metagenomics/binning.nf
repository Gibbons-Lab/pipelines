#!/usr/bin/env nextflow

params.data_dir = "data"

processed = Channel
    .fromFilePairs("${params.data_dir}/processed/*_R{1,2}_001.fastq.gz")
    .ifEmpty { error "Cannot find any read files in ${params.data_dir}!" }

assembled = Channel.from("{params.data_dir}/assembled/contigs/final.contigs.fa")


process assembly_align {
    cpus 16

    input:
    set id, file(reads) from processed
    file(assembly) from assembled

    output:
    set file("${id}.bam"), file("${id}.bai") into assembly_aligned

    """
    minimap2 -ax sr -t ${task.cpus} ${assembly} ${reads} | \
    samtools sort -@${task.cpus} -o ${id}.bam && \
    samtools index ${id}.bam ${id}.bai"
    """
}