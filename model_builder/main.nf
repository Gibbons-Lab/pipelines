#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${baseDir}/data"

params.media_db = null
params.media = null

max_threads = 24

def helpMessage() {
    log.info"""
    ~~~ Gibbons Lab Metabolic Model Builder Workflow ~~~

    Usage:
    A run using all,default parameters can be started with:
    > nextflow run main.nf --resume

    A run with all parametrs set would look like:
    > nextflow run main.nf --data_dir=./data --media_db=media.tsv --media="LB,M9"

    General options:
      --data_dir [str]              The main data directory for the analysis (must contain `raw`).
    Growth Media:
      --media_db                    A file containing growth media specification for CARVEME.
      --media                       Comma-separated list of media names to use.
    """.stripIndent()
}

params.help = false
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

process init_db {
  cpus 1
  publishDir "${params.data_dir}"

  output:
  path("db_stats.txt")

  script:
  """
  #! /usr/bin/env python

  import subprocess
  from carveme import config, project_dir
  from carveme.carve import first_run_check

  diamond_db = project_dir + config.get('generated', 'diamond_db')

  if name == "__main__":
    first_run_check()
    subprocess.call("diamond info -d %s > db_stats.txt" % diamond_db)
  """
}

process find_genes {
  cpus 1
  publishDir "${params.data_dir}/genes"

  input:
  tuple val(id), path(assembly)

  output:
  tuple val("${id}"), path("genes.ffn"), path("genes.faa")

  """
  prodigal -p meta -i ${assembly} -o genes.gff -d genes.ffn \
           -a genes.faa
  """
}

process build_model {
  cpus 1
  publishDir "${params.data_dir}/models"

  input:
  tuple val(id), path(genes_dna), path(genes_aa), path(db_info)

  output:
  tuple val("${id}"), path("${id}.xml.gz")

  script:
  if (params.media_db) && (params.media)
    """
    carve ${genes_aa} -o "${id}".xml.gz --mediadb ${params.media_db} --gapfill ${params.media}
    """
  else if (params.media)
    """
    carve ${genes_aa} -o "${id}".xml.gz --gapfill ${params.media}
    """
  else
    """
    carve ${genes_aa} -o "${id}".xml.gz
    """
}

process check_model {
  cpus 1
  publishDir "${params.data_dir}/model_qualities"

  input:
  tuple val(id), path(model)

  output:
  tuple val("${id}"), path("${id}.json")

  """
  memote run ${model} --ignore-git --location memote --filename ${id}.json
  """
}

workflow {
  Channel
    .fromPath([
        "${params.data_dir}/raw/*.fna",
        "${params.data_dir}/raw/*.fasta"
    ])
    .map{row -> tuple(row.baseName.split("\\.f")[0], path(row))}
    .set{genomes}

  init_db()
  find_genes(genomes)
  build_model(find_genes.out.combine(init_db.out))
  check_model(build_model.out)
}
