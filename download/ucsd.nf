nextflow.enable.dsl=2

params.barcodes = "barcode_map.csv"
params.url = "ftp://igm-storage.ucsd.edu"
params.run = ""
params.user = ""
params.password = ""


process get_run_files {
    publishDir "${launchDir}/data", mode: "copy", overwrite: true
    cpus 1

    output:
    path("listing.txt")

    """
    lftp ftp://${params.user}:${params.password}@igm-storage.ucsd.edu \
        -e "ls ${params.run}/; quit" > listing.txt
    """
}

process find_files {
    publishDir "${launchDir}/data", mode: "copy", overwrite: true
    cpus 1

    input:
    path(listing)
    path(barcodes)

    output:
    path("manifest.csv")

    """
    #!/usr/bin/env python

    import pandas as pd

    files = pd.read_csv("${listing}", sep=r"[\\t\\s]+", engine="python", header=None)
    files = files.iloc[:, [4, 8]]
    files.columns = ["size", "filename"]
    files["sample_name"] = files.filename.str.split("_S").str[0]

    barcodes = pd.read_csv("${barcodes}")
    barcodes.sample_name = barcodes.sample_name.str.replace(".", "_", regex=False)
    merged = barcodes.merge(files, on="sample_name")
    assert merged.shape[0] == 2 * barcodes.shape[0]  # paired-end

    merged.to_csv("manifest.csv", index=False)
    """
}

process download {
    publishDir "${launchDir}/data/raw", mode: "copy", overwrite: true
    cpus 24

    input:
    tuple val(id), val(filename)

    output:
    tuple val(id), path("*.fastq.gz")

    """
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 \
        ftp://${params.user}:${params.password}@igm-storage.ucsd.edu/${params.run}/${filename}
    """
}

workflow {
    Channel.fromPath("${launchDir}/data/${params.barcodes}").set{barcodes}

    get_run_files()

    find_files(get_run_files.out, barcodes)
    find_files.out
        .splitCsv(header : true)
        .map{row -> tuple(row.sample_name, row.filename)}
        .set{files}
    download(files)
}
