#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${baseDir}/data"
params.refs = "${baseDir}/refs"
params.eggnog_refs = "${params.refs}/eggnog"
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
params.threads = 12
params.ranks = "D,P,G,S"

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

Channel
    .of(params.ranks.split(","))
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
            --memory-mapping --report ${id}.tsv  ${reads[0]} ${reads[1]}
        """
}

process count_taxa {
    cpus 4

    input:
    tuple val(id), path(kraken), path(report), val(lev)

    output:
    tuple val(id), val(lev), path("${lev}/${id}.b2"), path("${lev}/${id}_bracken_mpa.tsv")

    """
    mkdir ${lev} && \
        sed 's/\\tR1\\t/\\tD\\t/g' ${report} > ${lev}/${report} && \
        bracken -d ${params.kraken2_db} -i ${lev}/${report} \
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
    from os import path
    import pandas as pd
    import re

    ranks = pd.Series({
        "d": "kingdom",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species"
    })

    def str_to_taxa(taxon):
        taxon = taxon.split("|")
        taxa = pd.Series({ranks[t.split("_", 1)[0]]: t.split("_", 1)[1] for t in taxon})
        return taxa

    read = []
    lev = "${lev}"
    input = "${reports.join(",")}"
    paths = input.split(",")

    for p in paths:
        print(f"processing {p}...")
        id = re.findall("(.+)_bracken_mpa.tsv", p)[0]
        try:
            counts = pd.read_csv(p, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            continue
        counts = counts[counts.iloc[:, 0].str.contains(
            str("d" if lev == "D" else lev).lower() + "_")]
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

process find_genes {
    cpus 1
    publishDir "${params.data_dir}/genes"

    input:
    tuple val(id), path(assembly)

    output:
    tuple path("${id}.ffn"), path("${id}.faa")

    """
    if grep -q ">" ${assembly}; then
        prodigal -p meta -i ${assembly} -o ${id}.gff -d ${id}.ffn -a ${id}.faa
    else
        touch ${id}.faa
        touch ${id}.ffn
    fi
    """
}

process cluster_transcripts {
    cpus params.threads/2
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(transcripts)

    output:
    path("transcripts.fna")

    """
    mmseqs easy-linclust ${transcripts} transcripts tmp \
        --cov-mode 0 -c ${params.overlap} \
        --min-seq-id ${params.identity} \
        --split-memory-limit 64G --threads ${task.cpus}
    rm transcripts_all_seqs.fna
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
    #!/usr/bin/env python

    from Bio import SeqIO
    import os

    os.system("cat ${proteins} > merged.faa")
    print("Reading transcript indices...")
    transcripts_idx = set(SeqIO.index("${transcripts}", "fasta"))
    print("Reading protein indices...")
    proteins = SeqIO.index("merged.faa", "fasta")
    print("Writing filtered proteins...")
    with open("proteins.faa", "wb") as out:
        for i, id in enumerate(transcripts_idx, start=1):
            out.write(proteins.get_raw(id))
            if (i % 100000) == 0:
                print(f"Processed {i} proteins.")
    os.system("rm merged.faa")
    """
}

process transcript_index {
    cpus params.threads
    publishDir "${params.data_dir}/"

    input:
    path(txns)

    output:
    path("txn_salmon_index")

    """
    salmon index -p ${task.cpus} -t ${txns} -i txn_salmon_index
    """
}

process map_and_count {
    cpus 8

    input:
    tuple val(id), path(reads), path(json), path(html), path(index)

    output:
    path("${id}.sf")

    script:
    if (params.single_end)
        """
        salmon quant --meta -p ${task.cpus} -l A -i ${index} -r ${reads} -o ${id} &&
            mv ${id}/quant.sf ${id}.sf || touch ${id}.sf
        """
    else
        """
        salmon quant --meta -p ${task.cpus} -l A -i ${index} -1 ${reads[0]} -2 ${reads[1]} -o ${id} &&
            mv ${id}/quant.sf ${id}.sf || touch ${id}.sf
        """
}

process merge_counts {
    cpus 1
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(salmon_quants)

    output:
    path("function_counts.csv.gz")

    """
    #!/usr/bin/env python

    from sys import stdin
    from os import path
    import pandas as pd
    import gzip

    paths = "${salmon_quants}"
    paths = paths.split(" ")
    with gzip.open("function_counts.csv.gz", "ab") as gzf:
        for i, p in enumerate(paths):
            sample = path.splitext(path.basename(p))[0]
            print("Processing sample {sample}...")
            try:
                counts = pd.read_csv(p, sep="\t").query("NumReads > 0.1")
            except Exception:
                continue

            counts.columns = [
                "locus_tag", "length", "effective_length", "tpm", "reads"]
            counts["sample_id"] = sample
            print(f"writing compressed output for sample {sample}...")
            counts.to_csv(gzf, header=(i==0),
                          index=False)
    """
}

process annotate {
    cpus params.threads
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(proteins)

    output:
    path("proteins.emapper.annotations")

    """
    TMP=\$(mktemp -d -t eggnog_results_XXXXXXXXXX)
    emapper.py -i ${proteins} --output proteins -m diamond \
        --data_dir ${params.eggnog_refs} --scratch_dir \$TMP --temp_dir /tmp \
        --cpu ${task.cpus}
    rm -rf \$TMP
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
                "${params.data_dir}/raw/*_{1,2}.fastq.gz",
                "${params.data_dir}/raw/*_R{1,2}.fastq.gz"
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
    //cluster_contigs(megahit.out.map{sample -> sample[1]}.collect())

    // find ORFs and collapse on 99% ANI
    find_genes(megahit.out)
    find_genes.out.map{sample -> sample[0]}.collect().set{transcripts}
    find_genes.out.map{sample -> sample[1]}.collect().set{proteins}
    cluster_transcripts(transcripts)
    filter_proteins(cluster_transcripts.out, proteins)

    // count gene abundances and annotate the genes
    transcript_index(cluster_transcripts.out)
    preprocess.out.combine(transcript_index.out) | map_and_count
    merge_counts(map_and_count.out.collect())
    annotate(filter_proteins.out)
}
