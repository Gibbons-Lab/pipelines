params.trim_left = 0
params.trunc_forward = 150
params.trunc_reverse = 150
params.maxEE = 2
params.merge = true
params.min_overlap = 8
params.forward_only = false
params.data_dir = "${baseDir}/data"
max_cpus = 24

process quality_control {
    publishDir "${params.data_dir}"
    cpus max_cpus

    output:
    tuple path("manifest.csv"), path("qc.rds"), path("qualities.png") into qc

    """
    #!/usr/bin/env Rscript
    library(mbtools)

    files <- find_read_files("${params.data_dir}/raw", dirs_are_runs = T)

    if (${params.forward_only ? "T" : "F"}) {
        files[, "reverse" := NULL]
    }

    fwrite(files, "manifest.csv")

    qc <- quality_control(files, min_score = 20)
    saveRDS(qc, "qc.rds")
    ggsave("qualities.png", pl = qc[["quality_plot"]] + theme_minimal(),
           width = 8, height = 4, dpi = 300)
    """
}

process trim {
    publishDir "${params.data_dir}"
    cpus max_cpus

    input:
    tuple path(manifest), path(qc), path(pl) from qc

    output:
    tuple path("preprocessed"), path("preprocessed.rds") into processed

    """
    #!/usr/bin/env Rscript
    library(mbtools)

    qc <- readRDS("${qc}")
    procced <- preprocess(
        qc,
        trimLeft = ${params.trim_left},
        truncLen = ${params.forward_only ? "${params.trunc_forward}" :
                     "c(${params.trunc_forward}, ${params.trunc_reverse})"},
        maxEE = ${params.maxEE},
        out_dir = "preprocessed",
        threads = ${task.cpus}
    )
    saveRDS(procced, "preprocessed.rds")
    """
}

process denoise {
    publishDir "${params.data_dir}"
    cpus max_cpus

    input:
    tuple path(procced), path(artifact) from processed

    output:
    tuple path("phyloseq.rds"), path("read_stats.csv"),
          path("denoised.rds") into denoised

    """
    #!/usr/bin/env Rscript
    library(mbtools)

    procced <- readRDS("${artifact}")
    procced[["files"]] <- procced[["files"]][procced[["passed"]][["preprocessed"]] > 0]
    denoised <- denoise(
        procced,
        hash = T,
        threads = ${task.cpus},
        merge = ${params.merge ? "T" : "F"},
        min_overlap = ${params.min_overlap}
    )
    saveRDS(denoised, "denoised.rds")
    fwrite(denoised[["passed_reads"]], "read_stats.csv")
    sdata <- as.data.frame(procced[["files"]])
    rownames(sdata) <- sdata[, "id"]
    ps <- as_phyloseq(denoised, sdata)
    saveRDS(ps, "phyloseq.rds")
    """
}

process tables {
    publishDir "${params.data_dir}"
    cpus 1

    input:
    tuple path(ps), path(stats), path(arti) from denoised

    output:
    tuple path("asvs.csv"), path("taxonomy.csv") into tabled

    """
    #!/usr/bin/env Rscript
    library(mbtools)

    denoised <- readRDS("${arti}")
    ids <- rownames(denoised[["feature_table"]])
    asvs <- as.data.table(denoised[["feature_table"]])[, "id" := ids]
    asvs <- melt(asvs, id.vars="id", variable.name="hash",
                 value.name="count")[count > 0]
    fwrite(asvs, "asvs.csv")
    ids <- rownames(denoised[["taxonomy"]])
    tax <- as.data.table(denoised[["taxonomy"]])[, "id" := ids]
    fwrite(tax, "taxonomy.csv")
    """
}
