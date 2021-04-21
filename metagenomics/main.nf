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
params.contig_length = 500
params.overlap = 0.8
params.identity = 0.99

max_threads = 24

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
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

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
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(taxonomy)

    output:
    path("multiqc_report.html")

    """
    multiqc ${params.data_dir}/preprocessed ${params.data_dir}/kraken2
    """
}


process megahit {
    cpus 4
    publishDir "${params.data_dir}/assembled"

    input:
    tuple val(id), path(reads), path(json), path(report)

    output:
    tuple val(id), path("contigs/${id}.contigs.fa")

    script:
    if (params.single_end)
        """
        megahit -r ${reads} -o contigs -t ${task.cpus} -m 0.4 \
                --min-contig-len ${params.contig_length} --out-prefix ${id}
        sed -i -e "s/^>/>${id}_/" contigs/${id}.contigs.fa
        """
    else
        """
        megahit -1 ${reads[0]} -2 ${reads[1]} -o contigs -t ${task.cpus} -m 0.4 \
                --min-contig-len ${params.contig_length} --out-prefix ${id}
        sed -i -e "s/^>/>${id}_/" contigs/${id}.contigs.fa
        """
}

process cluster_contigs {
    cpus max_threads/2
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(assemblies)

    output:
    path("all_contigs.fna")

    """
    cat ${assemblies} > merged.fna
    mmseqs easy-linclust merged.fna contigs tmp --cov-mode 0 -c ${params.overlap} --min-seq-id ${params.identity} --threads ${task.cpus}
    mv contigs_rep_seq.fasta all_contigs.fna
    """
}

process find_genes {
    cpus 1
    publishDir "${params.data_dir}/genes"

    input:
    tuple val(id), path(assembly)

    output:
    tuple path("${id}.ffn"), path("${id}.faa")

    """
    if grep -q ">" ${assembly}; then
        prodigal -p meta -i ${assembly} -o ${id}.gff -d ${id}.ffn \
             -a ${id}.faa
        sed -i -e "s/^>/>${id}_/" ${id}.faa
        sed -i -e "s/^>/>${id}_/" ${id}.ffn
    else
        touch ${id}.faa
        touch ${id}.ffn
    fi
    """
}

process cluster_transcripts {
    cpus max_threads/2
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(transcripts)

    output:
    path("transcripts.fna")

    """
    cat ${transcripts} > merged.fna
    mmseqs easy-linclust merged.fna transcripts tmp --cov-mode 0 -c ${params.overlap} --min-seq-id ${params.identity} --threads ${task.cpus}
    mv transcripts_rep_seq.fasta transcripts.fna
    """
}

process filter_proteins {
    cpus 1
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(transcripts)
    path(proteins)

    output:
    path("proteins.faa")

    """
    #!/usr/bin/env Rscript

    library(Biostrings)

    proteins <- "${proteins}"
    system2("cat", c(strsplit(proteins, " ")[[1]], ">", "merged.faa"))

    txns <- fasta.index("${transcripts}")
    txn_ids <- trimws(txns[["desc"]])
    prots <- readAAStringSet("merged.faa")
    prot_ids <- trimws(names(prots))
    names(prots) <- prot_ids
    writeXStringSet(prots[prot_ids %in% txn_ids], "proteins.faa")
    """
}

process transcript_index {
    cpus max_threads
    publishDir "${params.data_dir}/txn_bowtie2_index"

    input:
    path(txns)

    output:
    path("index")

    """
    mkdir index
    bowtie2-build -f --threads ${task.cpus} ${txns} index/txns
    """
}

process transcript_align {
    cpus 4
    publishDir "${params.data_dir}/txn_aligned"

    input:
    tuple val(id), path(reads), path(json), path(html), path(index)

    output:
    tuple val(id), path("${id}.bam")

    script:
    if (params.single_end)
        """
        bowtie2 -k 10 -p ${task.cpus} --mm --no-unal \
            -x ${index}/txns -1 ${reads[0]} -2 ${reads[1]} | \
            samtools view -bS - -o ${id}.bam
        """
    else
        """
        bowtie2 -k 10 -p ${task.cpus} --mm --no-unal \
            -x ${index}/txns -U ${reads} | \
            samtools view -bS - -o ${id}.bam
        """
}

process em_count {
    cpus 8

    input:
    tuple val(id), path(bam), path(genes)

    output:
    path("${id}.sf")

    """
    salmon quant -p ${task.cpus} -l SF -t ${genes} -a ${bam} -o ${id} &&
        mv ${id}/quant.sf ${id}.sf
    """
}

process merge_counts {
    cpus 1
    storeDir "${params.data_dir}", mode: "copy", overwrite: true

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
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(proteins)

    output:
    path("proteins.emapper.annotations")

    """
    emapper.py -i ${proteins} --output proteins -m diamond \
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
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

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
    // find files
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

    // quality filtering
    preprocess(raw)

    // quantify taxa abundances
    kraken(preprocess.out)
    count_taxa(kraken.out.combine(levels))
    count_taxa.out.map{s -> tuple(s[1], s[3])}
        .groupTuple()
        .set{merge_groups}
    merge_taxonomy(merge_groups)

    // quality overview
    multiqc(merge_taxonomy.out.collect())

    // assemble de novo
    megahit(preprocess.out)
    cluster_contigs(megahit.out.map{sample -> sample[1]}.collect())

    // find ORFs and collapse on 99% ANI
    find_genes(megahit.out)
    find_genes.out.map{sample -> sample[0]}.collect().set{transcripts}
    find_genes.out.map{sample -> sample[1]}.collect().set{proteins}
    cluster_transcripts(transcripts)
    filter_proteins(cluster_transcripts.out, proteins)

    // count gene abundances and annotate the genes
    transcript_index(cluster_transcripts.out)
    preprocess.out.combine(transcript_index.out) | transcript_align
    em_count(transcript_align.out.combine(cluster_transcripts.out))
    merge_counts(em_count.out.collect())
    annotate(filter_proteins.out)

     // contig_align(preprocess.out.combine(megahit.out))
    // replication_rates(contig_align.out)
    // merge_rates(replication_rates.out.collect())

}
