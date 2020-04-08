#!/usr/bin/env nextflow

params.data_dir = "data"
params.refs = "../refs/"
params.manifest = "manifest.csv"

params.trim_front = 5
params.min_length = 15
params.quality_threshold = 20
params.read_length = 100
params.threshold = 10

max_threads = 20

Channel
    .fromPath("${params.data_dir}/raw/*.fastq.gz")
    .ifEmpty { error "Cannot find any read files in ${params.data_dir}!" }
    .map{file -> tuple(file.baseName, file)}
    .set{raw}

Channel
    .from(["D", "P", "G", "S"])
    .set{levels, merge_levels}

process preprocess {
    cpus 1
    publishDir "${params.data_dir}/preprocessed"
    input:
    set id, file(reads) from raw

    output:
    set id, file("${id}.fastq.gz"),
            file("${id}.json"),
            file("${id}.html") into processed_assembly, processed_align,
                                    processed_kraken

    """
    fastp -i ${reads} -o ${id}.fastq.gz \
        --json ${id}.json --html ${id}.html \
        --trim_front1 ${params.trim_front} -l ${params.min_length} \
        -3 -M ${params.quality_threshold} -r ${reads}
    """
}

process kraken {
    cpus 4

    input:
    set id, file(forward), file(json), file(html) from processed_kraken

    output:
    set id, file("${id}.k2"), file("${id}.tsv") into kraken_reports

    """
    kraken2 --db /proj/gibbons/refs/kraken2_nt
        --threads ${task.cpus} --gzip-compressed --output ${id}.k2
        --report ${id}.tsv ${forward}"
    """
}

process count_taxa {
    cpus 4

    input:
    set id, file(kraken), file(report), lev from kraken_reports.combine(levels)

    output:
    set id, lev, file("${lev}/${id}.b2"), file("${lev}/${id}_mpa.tsv") into bracken_reports

    """
    cp ${report} ${lev}/${report} &&
        bracken -d /proj/gibbons/refs/kraken2_nt -i ${lev}/${report}
        -l ${lev} -o ${lev}/${id}.b2 -r {params.read_length}
        -t {params.threshold} &&
        kreport2mpa.py -r ${lev}/${id}_bracken.tsv -o ${lev}/${id}_mpa.tsv
        --no-intermediate-ranks
    """
}

bracken_reports
    .map{s -> tuple(s[1], s[3])}
    .groupTuple()
    .set{merge_groups}

process merge_taxonomy {
    cpus 1
    publishDir "${params.data_dir}"

    input:
    set lev, stdin file(reports)

    output:
    file("${lev}_counts.csv")

    """
    #!/usr/bin/env python

    from sys import stdin
    from loguru import logger as loggy
    from os import path
    import pandas as pd
    import re

    def str_to_taxa(taxon):
        taxon = taxon.split("|")
        taxa = pd.Series({t.split("_")[0]: t.split("_")[1] for t in taxon})
        return taxa

    read = []
    paths = stdin.readline().split(" ")

    for p in paths:
        lev, id = re.findall(r"/(\\w)/(\\w+)_mpa.tsv", p)
        counts = pd.read_csv(p, sep="\t", header=None)
        counts = counts[counts.iloc[:, 0].str.contains(
            str(lev).lower() + "_")]
        taxa = counts.iloc[:, 0].apply(str_to_taxa)
        taxa["reads"] = counts.iloc[:, 1]
        taxa["sample"] = id
        read.append(taxa)
    pd.concat(read, sort=False).to_csv("%s_counts.csv" % lev, index=False)
    """
}

merged = processed_assembly
    .map{it[1]}
    .reduce{a, b -> return "${a},${b}"}

process megahit {
    cpus max_threads
    publishDir "${params.data_dir}/assembled"

    input:
    val(forward) from merged

    output:
    file("contigs/final.contigs.fa") into (assembled_align, assembled_genes)

    """
    rm -rf contigs && \
    megahit -r ${forward} -o contigs -t ${task.cpus}
    """
}

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
    set id, file(forward), file(json), file(html), file(genes), file(aas) from processed_align.combine(genes)

    output:
    set id, file("${id}.bam") into transcript_aligned

    """
    minimap2 -acx sr -t ${task.cpus} ${genes} ${forward} | \
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
    conda "eggnog-mapper>=2.0.0"

    input:
    set file(genes), file(proteins) from genes_annotate

    output:
    file("denovo.emapper.annotations")

    """
    emapper.py -i ${proteins} --output denovo -m diamond \
        --data_dir ${params.refs}/eggnog-mapper/data \
        --cpu ${task.cpus} --resume
    """
}

