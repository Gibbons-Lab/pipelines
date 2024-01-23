nextflow.enable.dsl=2

params.runtable = "${baseDir}/data/runtable.csv"

process download {
    cpus 6
    maxRetries 3
    errorStrategy { (task.attempt <= maxRetries)  ? 'retry' : 'ignore' }
    scratch "/tmp"
    publishDir "${baseDir}/data/raw", mode: 'copy'

    input:
    val(run)

    output:
    path("*.fastq.gz")

    """
    aws s3 sync s3://sra-pub-run-odp/sra/${run} sra_${run} --no-sign-request
    fasterq-dump -e ${task.cpus} -f -3 ./sra_${run}/${run} && rm -rf ./sra_${run}
    pigz -p ${task.cpus} *.fastq
    """
}

process merge {
    publishDir "${baseDir}/data/merged"

    cpus 4

    input:
    tuple val(id), path(runs)

    output:
    path("${id}.fastq.gz")

    """
    cat ${runs} > ${id}.fastq.gz
    """
}

workflow {
    Channel.fromPath(params.runtable)
        .splitCsv(header: true)
        .map{row -> row.Run}
        .set{runs}

    runs | download
}
