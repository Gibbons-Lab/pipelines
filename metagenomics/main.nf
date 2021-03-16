#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${baseDir}/data"
params.refs = "/proj/gibbons/refs/"
params.eggnog_refs = "/tmp/eggnog"
params.kraken2_db = "${params.refs}/kraken2_default"

params.single_end = false
params.trim_front = 5
params.min_length = 50
params.quality_threshold = 20
params.read_length = 150
params.threshold = 10

max_threads = 24

Channel
    .fromList(["D", "P", "G", "S"])
    .set{levels}

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

process kraken {
    cpus 8
    publishDir "${params.data_dir}/kraken2"

    input:
    tuple val(id), path(reads), path(json), path(html)

    output:
    tuple val(id), path("${id}.k2"), path("${id}.tsv")

    script:
    if (params.single_end)
        """
        kraken2 --db ${params.kraken2_db} \
            --threads ${task.cpus} --gzip-compressed --output ${id}.k2 \
            --memory-mapping --report ${id}.tsv ${reads}
        """

    else
        """
        kraken2 --db ${params.kraken2_db} --paired \
            --threads ${task.cpus} --gzip-compressed --output ${id}.k2 \
            --memory-mapping --report ${id}.tsv ${reads[0]} ${reads[1]}
        """
}

process count_taxa {
    cpus 4

    input:
    tuple val(id), path(kraken), path(report), val(lev)

    output:
    tuple val(id), val(lev), path("${lev}/${id}.b2"), path("${lev}/${id}_bracken_mpa.tsv")

    """
    mkdir ${lev} && cp ${report} ${lev}/${report} && \
        bracken -d /proj/gibbons/refs/kraken2_default -i ${lev}/${report} \
        -l ${lev} -o ${lev}/${id}.b2 -r ${params.read_length} \
        -t ${params.threshold} -w ${lev}/${id}_bracken.tsv && \
        kreport2mpa.py -r ${lev}/${id}_bracken.tsv -o ${lev}/${id}_bracken_mpa.tsv \
        --no-intermediate-ranks
    """
}

process merge_taxonomy {
    cpus 1
    publishDir "${params.data_dir}"

    input:
    tuple val(lev), path(reports)

    output:
    path("${lev}_counts.csv")

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
    lev = "${lev}"
    input = "${reports.join(",")}"
    paths = input.split(",")

    for p in paths:
        loggy.warning(p)
        id = re.findall("(.+)_bracken_mpa.tsv", p)[0]
        try:
            counts = pd.read_csv(p, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            continue
        counts = counts[counts.iloc[:, 0].str.contains(
            str(lev).lower() + "_")]
        taxa = counts.iloc[:, 0].apply(str_to_taxa)
        taxa["reads"] = counts.iloc[:, 1]
        taxa["sample"] = id
        read.append(taxa.dropna())
    pd.concat(read, sort=False).to_csv("%s_counts.csv" % lev, index=False)
    """
}

process multiqc {
    publishDir "${params.data_dir}"

    input:
    path(taxonomy)

    output:
    path("multiqc_report.html")

    """
    multiqc ${params.data_dir}/preprocessed ${params.data_dir}/kraken2
    """
}


process megahit {
    cpus max_threads
    publishDir "${params.data_dir}/assembled"

    input:
    path(forward)
    path(reverse)

    output:
    path("contigs/final.contigs.fa")

    script:
    if (params.single_end)
        """
        megahit -r ${forward.join(",")} -o contigs -t ${task.cpus} -m 0.6
        """
    else
        """
        megahit -1 ${forward.join(",")} -2 ${reverse.join(",")} -o contigs -t ${task.cpus} -m 0.6
        """
}

process find_genes {
    cpus 1
    publishDir "${params.data_dir}/genes"

    input:
    path(assembly)

    output:
    tuple path("genes.ffn"), path("genes.faa")

    """
    prodigal -p meta -i ${assembly} -o genes.gff -d genes.ffn \
             -a genes.faa
    """
}

process transcript_align {
    cpus 4
    publishDir "${params.data_dir}/txn_aligned"

    input:
    tuple val(id), path(reads), path(json), path(html), path(genes), path(aas)

    output:
    tuple val(id), path("${id}.bam")

    """
    minimap2 -acx sr -t ${task.cpus} ${genes} ${reads} | \
    samtools view -bS - -o ${id}.bam
    """
}

process em_count {
    cpus 8

    input:
    tuple val(id), path(bam), path(genes), path(aas)

    output:
    path("${id}.sf")

    """
    salmon quant -p ${task.cpus} -l SF -t ${genes} -a ${bam} -o ${id} &&
        mv ${id}/quant.sf ${id}.sf
    """
}

process merge_counts {
    cpus 1
    storeDir "${params.data_dir}"

    input:
    path(salmon_quants)

    output:
    path("function_counts.csv.gz")

    """
    #!/usr/bin/env python

    from sys import stdin
    from loguru import logger as loggy
    from os import path
    import pandas as pd

    read = []
    paths = "${salmon_quants}"
    paths = paths.split(" ")
    for p in paths:
        sample = path.splitext(path.basename(p))[0]
        loggy.info("Processing sample {}...", sample)
        counts = pd.read_csv(p, sep="\t").query("NumReads > 0.1")

        counts.columns = [
            "locus_tag", "length", "effective_length", "tpm", "reads"]
        counts["sample"] = sample
        read.append(counts)
    read = pd.concat(read)
    loggy.info("writing compressed output")
    read.to_csv("function_counts.csv.gz", index=False)
    """
}

process annotate {
    cpus max_threads
    publishDir "${params.data_dir}/annotated"

    input:
    tuple path(genes), path(proteins)

    output:
    path("denovo.emapper.annotations")

    """
    set +eu && \
    source ~/miniconda3/etc/profile.d/conda.sh && \
    conda activate eggnog && \
    emapper.py -i ${proteins} --output denovo -m diamond \
        --data_dir ${params.eggnog_refs} \
        --cpu ${task.cpus} --resume
    """
}

process contig_align {
    cpus 4
    publishDir "${params.data_dir}/contig_aligned"

    input:
    tuple val(id), path(reads), path(json), path(html), path(contigs)

    output:
    tuple val(id), path("${id}.bam")

    """
    minimap2 -acx sr -t ${task.cpus} ${contigs} ${reads} | \
    samtools view -bS - -o ${id}.bam
    """
}

process replication_rates {
    cpus 8
    publishDir "${params.data_dir}/replication_rates"

    input:
    tuple val(id), path(bam)

    output:
    path("${id}.csv")

    """
    #!/usr/bin/env Rscript

    library(mbtools)
    library(futile.logger)

    flog.threshold(DEBUG)
    alns <- data.table(id=${id}, alignment=${bam}, success=TRUE)

    co <- bin_coverage(alns, threads=${task.cpus})
    rates <- replication_rates(co, threads=${task.cpus}, min_points_fit=40)
    fwrite(rates[["rate"]], "${id}.csv")
    """
}

process merge_rates {
    cpus 1
    publishDir "${params.data_dir}"

    input:
    path(rates)

    output:
    path("replication_rates.csv")

    """
    #!/usr/bin/env Rscript

    library(data.table)

    files <- strsplit("${rates}", " ")[[1]]
    rates <- rbindlist(lapply(files, fread))
    fwrite(files, "replication_rates.csv")
    """
}

workflow {
    if (params.single_end) {
        Channel
            .fromPath("${params.data_dir}/raw/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{raw}
    } else {
        Channel
            .fromFilePairs([
                "${params.data_dir}/raw/*_R{1,2}_001.fastq.gz",
                "${params.data_dir}/raw/*_{1,2}.fastq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${params.data_dir}!" }
            .set{raw}
    }

    preprocess(raw)

    kraken(preprocess.out)
    count_taxa(kraken.out.combine(levels))
    count_taxa.out.map{s -> tuple(s[1], s[3])}
        .groupTuple()
        .set{merge_groups}
    merge_taxonomy(merge_groups)

    multiqc(merge_taxonomy.out.collect())

    if (params.single_end) {
        forward = preprocess.out.map{it -> it[1]}.toSortedList()
        megahit(forward, forward)
    } else {
        forward = preprocess.out.map{it -> it[1][0]}.toSortedList()
        reverse = preprocess.out.map{it -> it[1][1]}.toSortedList()
        megahit(forward, reverse)
    }

    find_genes(megahit.out)
    preprocess.out.combine(find_genes.out) | transcript_align
    em_count(transcript_align.out.combine(find_genes.out))
    merge_counts(em_count.out.collect())

    contig_align(preprocess.out.combine(megahit.out))
    replication_rates(contig_align.out)
    merge_rates(replication_rates.out.collect())

    annotate(find_genes.out)
}
