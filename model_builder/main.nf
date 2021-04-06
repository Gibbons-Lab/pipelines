#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${baseDir}/data"
params.media_db = null
params.media = null


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
  publishDir "${params.data_dir}", mode: "copy", overwrite: true

  output:
  path("db_stats.txt")

  """
  #!/usr/bin/env python

  import subprocess
  from carveme import config, project_dir
  from carveme.cli.carve import first_run_check

  diamond_db = project_dir + config.get('generated', 'diamond_db')[:-5] + ".dmnd"

  if __name__ == "__main__":
    first_run_check()
    with open("db_stats.txt", "w") as out:
      res = subprocess.check_output(["diamond", "dbinfo", "-d", diamond_db])
      out.write(res.decode())
  """
}

process download_bigg {
  cpus 1

  output:
  path("universal_model.json")

  """
  wget http://bigg.ucsd.edu/static/namespace/universal_model.json
  """
}

process find_genes {
  cpus 1
  publishDir "${params.data_dir}/genes", mode: "copy", overwrite: true

  input:
  tuple val(id), path(assembly)

  output:
  tuple val("${id}"), path("${id}.ffn"), path("${id}.faa")

  """
  prodigal -p single -i ${assembly} -o ${id}.gff -d ${id}.ffn \
           -a ${id}.faa
  """
}

process build_model {
  cpus 2

  input:
  tuple val(id), path(genes_dna), path(genes_aa), path(db_info)

  output:
  tuple val("${id}"), path("${id}_draft.xml.gz")

  script:
  if (params.media_db && params.media)
    """
    carve ${genes_aa} -o "${id}_draft".xml.gz --mediadb ${params.media_db} \
          --gapfill ${params.media} --diamond-args "-p ${task.cpus}" --fbc2 -v
    """
  else if (params.media)
    """
    carve ${genes_aa} -o "${id}_draft".xml.gz --gapfill ${params.media} \
          --diamond-args "-p ${task.cpus}" --fbc2 -v
    """
  else
    """
    carve ${genes_aa} -o "${id}_draft".xml.gz --diamond-args "-p ${task.cpus}" --fbc2 -v
    """
}

process annotate_model {
  cpus 1
  publishDir "${params.data_dir}/models", mode: "copy", overwrite: true

  input:
  tuple val(id), path(model), path(bigg)

  output:
  tuple val("${id}"), path("${id}.xml.gz")

  """
  #!/usr/bin/env python

  import json
  import pandas as pd
  from collections import defaultdict
  from cobra.io import read_sbml_model, write_sbml_model

  def parse_annotation(anns):
    parsed = defaultdict(list)
    for _, url in anns:
      key, id = url.replace("http://identifiers.org/", "").split("/", 1)
      parsed[key] += [id]
    return dict(parsed)

  bigg = json.load(open("${bigg}"))
  reactions = pd.DataFrame.from_records(bigg["reactions"])
  reactions.index = reactions.id
  metabolites = pd.DataFrame.from_records(bigg["metabolites"])
  metabolites.index = metabolites.id

  model = read_sbml_model("${model}")

  for m in model.metabolites:
    m.formula = m.formula.split(";")[0]

  reactions = reactions[reactions.index.isin(r.id for r in model.reactions)]
  metabolites = metabolites[metabolites.index.isin(m.id for m in model.metabolites)]

  for rid, entries in reactions.iterrows():
    rxn = model.reactions.get_by_id(rid)
    rxn.annotation = parse_annotation(entries["annotation"])

  for mid, entries in metabolites.iterrows():
    met = model.metabolites.get_by_id(mid)
    met.annotation = parse_annotation(entries["annotation"])

  write_sbml_model(model, "${id}.xml.gz")
  """
}

process check_model {
  cpus 1
  publishDir "${params.data_dir}/model_qualities", mode: "copy", overwrite: true

  input:
  tuple val(id), path(model)

  output:
  tuple val("${id}"), path("${id}.html")

  """
  memote report snapshot ${model} --solver cplex --filename ${id}.html
  """
}

workflow {
  Channel
    .fromPath([
        "${params.data_dir}/raw/*.fna",
        "${params.data_dir}/raw/*.fasta"
    ])
    .map{row -> tuple(row.baseName.split("\\.f")[0], tuple(row))}
    .set{genomes}

  init_db()
  download_bigg()
  find_genes(genomes)
  build_model(find_genes.out.combine(init_db.out))
  annotate_model(build_model.out.combine(download_bigg.out))
  check_model(annotate_model.out)
}
