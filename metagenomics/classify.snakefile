import pandas as pd
from os import path
from tempfile import TemporaryDirectory

configfile: "config.yml"

manifest = pd.read_csv(path.join(config["data"], config["manifest"]))
samples = manifest.iloc[:, 1].str.split("_R\\d_001.+$").str[0]
tempdir = TemporaryDirectory("snakemake")

def str_to_taxa(taxon):
    taxon = taxon.split("|")
    taxa = pd.Series({t.split("_")[0]: t.split("_")[1] for t in taxon})
    return taxa

def summarize_rank(counts, level):
    lidx = counts.columns.get_loc(level)
    valid = counts.iloc[:, lidx].notnull()
    valid &= counts.iloc[:, (lidx+1):8].apply(
        lambda s: pd.isnull(s).all(),
        axis=1)
    summarized = (
        counts[valid].groupby(counts.columns[0:(lidx + 1)].
        tolist() + ["sample"]).reads.sum().reset_index()
    )
    return summarized

rule all:
    input:
        expand(
            path.join(config["data"], "{lev}_counts.csv"),
            lev=["D", "P", "G", "S"]
        )

rule classify:
    input:
        path.join(config["ref"], "kraken2_extended"),
        path.join(config["data"], "preprocessed", "{sample}_R1_001.fastq.gz"),
        path.join(config["data"], "preprocessed", "{sample}_R2_001.fastq.gz")
    output:
        path.join(config["data"], "taxonomy", "{sample}.k2"),
        path.join(config["data"], "taxonomy", "{sample}.tsv")
    threads: 20
    conda: "conda_classify.yml"
    shell:
        "kraken2 --db {input[0]} "
        "--threads {threads} --gzip-compressed --output {output[0]} "
        " --report {output[1]} --paired {input[1]} {input[2]}"

rule count:
    input:
        path.join(config["data"], "taxonomy", "{sample}.tsv"),
        path.join(config["ref"], "kraken2_extended")
    output:
        path.join(config["data"], "taxonomy", "{level}", "{sample}.b2"),
        path.join(config["data"], "taxonomy", "{level}", "{sample}_mpa.tsv")
    params:
        threshold = 10,
        read_length = 100,
        bracken_input = path.join(config["data"], "taxonomy", "{level}", "{sample}.tsv"),
        bracken_report = path.join(config["data"], "taxonomy", "{level}", "{sample}_bracken.tsv")
    threads: 1
    conda: "conda_classify.yml"
    shell:
        "cp -f {input[0]} {params.bracken_input} && "
        "bracken -d {input[1]} -i {params.bracken_input} "
        "-l {wildcards.level} -o {output[0]} -r {params.read_length} "
        "-t {params.threshold} && "
        "kreport2mpa.py -r {params.bracken_report} -o {output[1]} "
        "--no-intermediate-ranks"

rule merge:
    input:
        expand(
            path.join(config["data"], "taxonomy", "{{level}}", "{sample}.b2"),
            sample=samples
        )
    output:
        path.join(config["data"], "{level}_counts.csv")
    threads: 1
    run:
        taxa_files = []
        for sa in samples:
            counts = pd.read_csv(
                path.join(config["data"], "taxonomy", wildcards.level,
                          "%s_mpa.tsv" % (sa)),
                sep="\t",
                header=None
            )
            counts = counts[counts.iloc[:, 0].str.contains(
                str(wildcards.level).lower() + "_")]
            taxa = counts.iloc[:, 0].apply(str_to_taxa)
            taxa["reads"] = counts.iloc[:, 1]
            taxa["sample"] = sa
            taxa_files.append(taxa)
        pd.concat(taxa_files, sort=False).to_csv(str(output), index=False)

