nextflow.enable.dsl=2

params.runtable = "${baseDir}/data/runtable.csv"

process download {
    cpus 6
    publishDir "${baseDir}/data/raw", mode: 'copy'

    input:
    val(run)

    output:
    path("*.fastq.gz")

    """
    fasterq-dump -e ${task.cpus} -f -3 ${run}
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
