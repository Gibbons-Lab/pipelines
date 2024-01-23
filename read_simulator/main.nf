#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.read_length = 150
params.genomes = "${launchDir}/data/genomes"
params.mutation_rate = 0.0
params.manifest = "${launchDir}/data/manifest.csv"

process sample {
    cpus 1
    memory "4 GB"
    beforeScript "ulimit -Sf unlimited"

    input:
    tuple val(sample_id), val(id), path(genome), val(n)

    output:
    tuple val(sample_id), val(id), path("*.fastq.gz")

    """
    trap "rm ref.fna" EXIT

    zcat ${genome} > ref.fna

    dwgsim -N ${n} -H \
        -c 0 -S 0 -A 0 -e "0.001-0.005" -E "0.005-0.01" \
        -1 ${params.read_length} -2 ${params.read_length} -n ${params.read_length} \
        -d ${params.read_length * 2} \
        -r ${params.mutation_rate} -y 0 \
        -P ${id} -o 1 ref.fna sim

    mv sim.bwa.read1.fastq.gz ${sample_id}_${id}_R1.fastq.gz
    mv sim.bwa.read2.fastq.gz ${sample_id}_${id}_R2.fastq.gz
    """
}

process merge_reads {
    cpus 1
    memory "4 GB"

    input:
    tuple val(sample_id), path(forward), path(reverse)

    output:
    tuple val(sample_id), path("${sample_id}_*.fastq.gz")

    """
    cat ${forward} > ${sample_id}_R1.fastq.gz
    cat ${reverse} > ${sample_id}_R2.fastq.gz
    """

}

process random_order {
    cpus 4
    memory "16 GB"
    scratch "/tmp"
    publishDir "${launchDir}/data/sampled", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("shuffled/*.fastq.gz")

    """
    #!/usr/bin/env python

    from Bio import SeqIO
    import gzip
    import os
    from random import shuffle
    from pathlib import Path
    import shutil

    base = Path("shuffled")

    with gzip.open("${reads[0]}", "rt") as infile, open("fwd.fastq", "w") as outfile:
        shutil.copyfileobj(infile, outfile)
    with gzip.open("${reads[1]}", "rt") as infile, open("rev.fastq", "w") as outfile:
        shutil.copyfileobj(infile, outfile)
    fidx = SeqIO.index("fwd.fastq", "fastq")
    ridx = SeqIO.index("rev.fastq", "fastq")
    print(f"Finished indexing {len(fidx)} records...")
    ord = list(range(len(fidx)))
    shuffle(ord)

    os.mkdir(base)
    print("Sampling records...", flush=True)

    fnames = list(fidx.keys())
    rnames = list(ridx.keys())
    with gzip.open(base / "${sample_id}_R1.fastq.gz", "wb") as fwd, gzip.open(base / "${sample_id}_R2.fastq.gz", "wb") as rev:
        for i, rec in enumerate(ord):
            fwd.write(fidx.get_raw(fnames[rec]))
            rev.write(ridx.get_raw(rnames[rec]))
            if i % 10000 == 0:
                print(f"Processed {i} records...", flush=True)

    os.remove("fwd.fastq")
    os.remove("rev.fastq")
    print("Done.")
    """

}

def relative_depth(abundance, depth) {
    Math.round(Float.parseFloat(abundance) * Float.parseFloat(depth))
}

workflow {
    Channel
        .fromPath("${params.manifest}")
        .splitCsv(header: true)
        .map{row -> tuple(
                row.sample_id, row.id, "${params.genomes}/${row.file}",
                relative_depth(row.relative_abundance, row.depth)
            )
        }
        .set{manifest}

    manifest | sample
    sample.out
        .map{ it -> tuple(it[0], it[2][0]) }
        .groupTuple()
        .set{forward_groups}
    sample.out
        .map{ it -> tuple(it[0], it[2][1]) }
        .groupTuple()
        .set{reverse_groups}
    forward_groups
        .combine(reverse_groups, by: 0)
        .set{fastq_groups}

    fastq_groups | merge_reads | random_order

}
