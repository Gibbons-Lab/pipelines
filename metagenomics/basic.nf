#!/usr/bin/env nextflow

params.data_dir = "data"
params.refs = "../refs/"
params.manifest = "manifest.csv"

params.trim_front = 5
params.min_length = 15
params.quality_threshold = 20

max_threads = 32

raw = Channel
    .fromFilePairs("${params.data_dir}/raw/*_R{1,2}_001.fastq.gz")
    .ifEmpty { error "Cannot find any read files in ${params.data_dir}!" }

process preprocess {
    cpus 1
    publishDir "${params.data_dir}/preprocessed"
    input:
    set id, file(reads) from raw

    output:
    set id, file("${id}_R1_001.fastq.gz"),
            file("${id}_R2_001.fastq.gz"),
            file("${id}.json"),
            file("${id}.html") into processed_assembly, processed_align

    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${id}_R1_001.fastq.gz -O ${id}_R2_001.fastq.gz \
        --json ${id}.json --html ${id}.html \
        --trim_front1 ${params.trim_front} -l ${params.min_length} \
        -3 -M ${params.quality_threshold} -r ${reads[0]}
    """
}

process megahit {
    cpus max_threads
    publishDir "${params.data_dir}/assembled"

    input:
    val(all) from processed_assembly.collect()

    output:
    file("contigs/final.contigs.fa") into (assembled_align, assembled_genes)

    def reads = all
        .map{tuple(it[1], it[2])}
        .reduce{a, b -> return tuple("${a[0]},${b[0]}", "${a[1]},${b[1]}")}

    """
    rm -rf contigs && \
    megahit -1 ${all[0]} -2 ${all[1]} -o contigs -t ${task.cpus}
    """
}

Channel
    .from(file("${params.data_dir}/assembled/contigs/final.contigs.fa"))
    .into{assembled_align; assembled_genes}

process find_genes {
    cpus 1
    publishDir "${params.data_dir}/genes"

    input:
    file(assembly) from assembled_genes

    output:
    set file("genes.ffn"), file("genes.faa") into (genes, genes_count, genes_annotate)

    """
    prodigal -p meta -i ${assembly} -o genes.gff -d genes.ffn \
             -a genes.faa
    """
}

process transcript_align {
    cpus 8
    publishDir "${params.data_dir}/txn_aligned"

    input:
    set id, file(forward), file(reverse), file(json), file(html), file(genes), file(aas) from processed_align.combine(genes)

    output:
    set id, file("${id}.bam") into transcript_aligned

    """
    minimap2 -acx sr -t ${task.cpus} ${genes} ${forward} ${reverse} | \
    samtools view -bS - -o ${id}.bam
    """
}

process em_count {
    cpus max_threads

    input:
    set id, file(bam), file(genes), file(aas) from transcript_aligned.combine(genes_count)

    output:
    file("${id}/quant.sf") into salmon_quants

    """
    salmon quant -p ${task.cpus} -l SF -t ${genes} -a ${bam} -o ${id} &> /dev/null
    """
}

process merge {
    cpus 1
    storeDir "${params.data_dir}"

    input:
    stdin salmon_quants.reduce{a, b -> "$a,$b"}

    output:
    file("function_counts.csv.gz")

    """
    #!/usr/bin/env python

    from sys import stdin
    from loguru import logger as loggy
    from os import path
    import pandas as pd

    read = []
    paths = stdin.readline().split(",")
    for p in paths:
        sample = path.basename(path.dirname(p))
        loggy.info("Processing sample {}...", sample)
        counts = pd.read_csv(p, sep="\t").query("NumReads > 0.1")

        counts.columns = [
            "locus_tag", "length", "effective_length", "tpm", "reads"]
        counts["sample"] = sample
        read.append(counts)
    read = pd.concat(read)
    loggy.info("writing compressed output")
    read.to_csv("function_counts.csv.gz")
    """
}

process annotate {
    cpus max_threads
    publishDir "${params.data_dir}/annotated"

    input:
    set file(genes), file(proteins) from genes_annotate

    output:
    file("denovo.emapper.annotations")

    """
    conda activate eggnog && \
    emapper.py -i ${proteins} --output denovo -m diamond \
        --data_dir /proj/gibbons/refs/eggnog-mapper/data \
        --cpu ${task.cpus} --resume
    """
}
