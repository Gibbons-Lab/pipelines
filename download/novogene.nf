nextflow.enable.dsl=2

params.runtable = "${baseDir}/data/links.txt"

process download {
    cpus 6
    publishDir "${baseDir}/data/raw", mode: 'copy'

    input:
    val(link)

    output:
    path("*.fq.gz")

    """
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 ${link}
    """
}

workflow {
    Channel.fromPath(params.runtable)
        .splitCsv()
        .map{ row -> row[0] }
        .filter( s -> s=~/RawData.*\.fq\.gz$/ )
        .set{runs}

    runs | download
}
