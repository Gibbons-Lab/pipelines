params.trim_left = 0
params.trunc_forward = 240
params.trunc_reverse = 230
params.maxEE = 2
params.merge = true
params.min_overlap = 8
params.forward_only = false
params.data_dir = "${baseDir}/data"
params.taxa_db = "/proj/gibbons/refs/silva_nr99_v138.1_train_set.fa.gz"
params.species_db = "/proj/gibbons/refs/silva_species_assignment_v138.1.fa.gz"
params.threads = 16
params.pattern = "illumina"

process quality_control {
    publishDir "${params.data_dir}", mode: "copy", overwrite: true
    cpus params.threads

    output:
    tuple path("manifest.csv"), path("qc.rds"), path("*.png")

    """
    #!/usr/bin/env Rscript
    library(mbtools)

    files <- find_read_files(
        "${params.data_dir}/raw",
        pattern = mbtools:::${params.pattern}_pattern,
        annotations = mbtools:::${params.pattern}_annotations,
        dirs_are_runs = T
    )

    if (${params.forward_only ? "T" : "F"}) {
        files[, "reverse" := NULL]
    }

    fwrite(files, "manifest.csv")

    qc <- quality_control(files, min_score = 20)
    saveRDS(qc, "qc.rds")
    ggsave("qualities.png", pl = qc[["quality_plot"]] + theme_minimal(),
           width = 8, height = 4, dpi = 300)
    ggsave("entropy.png", pl = qc[["entropy_plot"]] + theme_minimal(),
           width = 8, height = 4, dpi = 300)
    ggsave("lengths.png", pl = qc[["length_plot"]] + theme_minimal(),
           width = 8, height = 4, dpi = 300)
    """
}

process trim {
    publishDir "${params.data_dir}"
    cpus params.threads

    input:
    tuple path(manifest), path(qc), path(pl)

    output:
    tuple path("preprocessed"), path("preprocessed.rds")

    """
    #!/usr/bin/env Rscript
    library(mbtools)

    qc <- readRDS("${qc}")
    manifest <- fread("${manifest}")

    if ("reverse" %in% names(manifest)) {
        trunc <- c(${params.trunc_forward}, ${params.trunc_reverse})
    } else {
        trunc <- ${params.trunc_forward}
    }

    procced <- preprocess(
        qc,
        trimLeft = ${params.trim_left},
        truncLen = trunc,
        maxEE = ${params.maxEE},
        out_dir = "preprocessed",
        threads = ${task.cpus}
    )
    saveRDS(procced, "preprocessed.rds")
    """
}

process denoise {
    publishDir "${params.data_dir}/denoise"
    cpus params.threads

    input:
    tuple path(procced), path(artifact)

    output:
    tuple path("read_stats.csv"), path("denoised.rds"), path("phyloseq.rds")

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
        min_overlap = ${params.min_overlap},
        taxa_db = "${params.taxa_db}",
        species_db = "${params.species_db}"
    )
    saveRDS(denoised, "denoised.rds")
    fwrite(denoised[["passed_reads"]], "read_stats.csv")
    sdata <- as.data.frame(procced[["files"]])
    rownames(sdata) <- sdata[, "id"]
    ps <- as_phyloseq(denoised, sdata)
    saveRDS(ps, "phyloseq.rds")
    """
}

process tree {
    publishDir "${params.data_dir}", mode: 'copy', overwrite: true

    cpus params.threads

    input:
    tuple path(stats), path(denoised), path(ps)

    output:
    tuple path("asvs.tree"), path("phyloseq.rds")

    """
    #!/usr/bin/env Rscript

    library(futile.logger)
    library(mbtools)

    denoised <- readRDS("${denoised}")

    seqs <- denoised[["taxonomy"]][, "sequence"]
    alignments <- DECIPHER::AlignSeqs(
        Biostrings::DNAStringSet(seqs),
        anchor = NA,
        processors = ${task.cpus}
    )

    am <- as(alignments, "matrix")
    freqs <- t(apply(am, 2, function(x)
        table(factor(x, levels = c("-", "A", "C", "G", "T"))) / length(x)))
    absent <- freqs[, "-"] > 0.5
    flog.info(paste("%d alignment columns are absent in >50%% of alignments,",
                    "removing them."), sum(absent))
    alignments <- am[, !absent] %>% apply(1, paste0, collapse = "") %>%
                  Biostrings::DNAStringSet()
    flog.info("Final alignment has %d positions. Starting FastTree...",
              width(alignments)[1])
    Biostrings::writeXStringSet(alignments, "asvs_aligned.fna")
    args <- c("-fastest", "-gtr", "-gamma",
              "-nt", "asvs_aligned.fna", ">", "asvs.tree")
    out <- system2("FastTreeMP",
                   args = args,
                   env = "OMP_NUM_THREADS=${task.cpus}"
            )
    tree <- read_tree("asvs.tree")

    ps <- readRDS("${ps}")
    phy_tree(ps) <- tree

    saveRDS(ps, "phyloseq.rds")
    """
}


process tables {
    publishDir "${params.data_dir}", mode: "copy", overwrite: true
    cpus 1

    input:
    tuple path(stats), path(arti), path(ps)

    output:
    tuple path("asvs.csv"), path("taxonomy.csv")

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

workflow {
    quality_control | trim | denoise | tables
    denoise.output | tree
}
