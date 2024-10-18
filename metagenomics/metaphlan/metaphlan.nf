#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "/proj/gibbons/data"
params.fastq_dir="${params.data_dir}/fastq"

params.mpa_db ="/proj/gibbons/refs/metaphlan4_db_202403"
params.mpa_index="mpa_vJun23_CHOCOPhlAnSGB_202403"

params.single_end = false
params.trim_front = 5
params.min_length = 50
params.quality_threshold = 20
params.read_length = 150
params.threads = 12

def helpMessage() {
    log.info"""
    ~~~ Reference Based Metagenomics Workflow with Metaphlan4 ~~~

    Usage:
    A run using all,default parameters can be started with:
    > nextflow run ref_analyis.nf -resume

    A run with all parametrs set would look like:
    > nextflow run metaphlan_analysis.nf --data_dir=./data/analysis --fastq_dir=./data/fastq --single_end=false \\
                           --trim_front=5 --min_length=50 --quality_threshold=20 --read_length=150 --threads=12

    General options:
      --data_dir [str]              The main data directory for the analysis.
      --fastq_dir [str]             Directory containing raw data.
      --read_length [str]           The length of the reads.
      --single_end [bool]           Specifies that the input is single-end reads.
      --threads [int]               The maximum number of threads a single process can use.
                                    This is not the same as the maximum number of total threads used.
    Reference DBs:
      --refs [str]                  Folder in which to find references DBs.

    Quality filter:
      --trim_front [str]            How many bases to trim from the 5' end of each read.
      --min_length [str]            Minimum accepted length for a read.
      --quality_threshold [str]     Smallest acceptable average quality.
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

process metaphlan {
    cpus 5
    publishDir "${params.data_dir}/metaphlan", mode:'copy', overwrite: true

    input:
    tuple val(id), path(reads), path(json), path(html)

    output:
    tuple val(id), file("${id}_metaphlan.txt")

    script:
    "metaphlan ${params.data_dir}/preprocessed/${reads[0]},${params.data_dir}/preprocessed/${reads[1]} --bowtie2out ${id}_metagenome.bowtie2.bz2 --nproc ${task.cpus} --input_type fastq -o ${id}_metaphlan.txt --index ${params.mpa_index} --bowtie2db ${params.mpa_db} --unclassified_estimation --add_viruses "
}

process merge_metaphlan {
    publishDir "${params.data_dir}/metaphlan", mode:'copy', overwrite: true
    
    input:
    val(profiles)

    output:
    file("metaphlan_abundance_table.txt")
    
    script:
    "merge_metaphlan_tables.py ${params.data_dir}/metaphlan/*_metaphlan.txt > metaphlan_abundance_table.txt"

}

process reformat_output {
    publishDir "${params.data_dir}/metaphlan", mode:'copy', overwrite: true

    input:
    file("metaphlan_abundance_table.txt")

    output:
    tuple file("P_counts.csv"), file("F_counts.csv"), file("G_counts.csv"), file("S_counts.csv")

    """
    #!/usr/bin/env python

    import pandas as pd
    import numpy as np
    
    def get_tax(df,level='g',taxcol='clade_name'):
        tax='kpcofgst'
        pos = tax.index(level)
        if pos < len(tax):
            df=df[(df[taxcol].str.contains('%s__'%(level)))&(~df.clade_name.str.contains('%s__'%(tax[pos+1])))]
        else:
            df=df[df[taxcol].str.contains('%s__'%(level))]
        for i in range(pos+1):
            tx=tax[i]
            df.loc[:,tx]=np.nan
        for i in df.index:
            for tx in df.loc[i,taxcol].split('|'):
                df.loc[i,tx[0]]=tx.split('__')[1]
        df=df.drop(taxcol,axis=1)
        return df

    def reformat(df):
        samples=df.columns[df.dtypes!=object]
        tax=df.columns[df.dtypes==object]
        tax_dict=df[tax].to_dict()
        df=df[samples].stack().reset_index()
        df.columns=['keys','sample','reads']
        for col in tax_dict.keys():
            d=tax_dict[col]
            df[col]=df['keys'].map(d)
        df=df.drop('keys',axis=1)
        return df.sort_values(by=['reads'],ascending=False)

    metaphlan=pd.read_csv('metaphlan_abundance_table.txt',sep='\t',skiprows=1)
    metaphlan.columns = [x if 'metaphlan' not in x else x[:-10] for x in metaphlan.columns]
    for tax in 'PFGS':
        tx = tax.lower()
        df=get_tax(metaphlan,level=tx)
        df=reformat(df)
        df.to_csv('%s_counts.csv'%tax)

    """

}

process multiqc {
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    file("metaphlan_abundance_table.txt")

    output:
    path("multiqc_report.html")

    """
    multiqc ${params.data_dir}/preprocessed ${params.data_dir}/metaphlan
    """
}

workflow {
    // find files
    if (params.single_end) {
        Channel
            .fromPath("${params.fastq_dir}/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{raw}
    } else {
        Channel
            .fromFilePairs([
                "${params.fastq_dir}/*_R{1,2}_001.fastq.gz",
                "${params.fastq_dir}/*_{1,2}.fastq.gz",
                "${params.fastq_dir}/*_R{1,2}.fastq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${params.fastq_dir}!" }
            .set{raw}
    }

    // quality filtering
    preprocess(raw)
    // quantify taxa abundances
    metaphlan(preprocess.out)
    // merge results
    merge_metaphlan(metaphlan.out.collect())
    // reformat and extract data for specific taxonomic levels
    reformat_output(merge_metaphlan.out)
    // generate html report
    multiqc(merge_metaphlan.out)
}
