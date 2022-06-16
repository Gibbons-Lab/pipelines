#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${baseDir}/data"
params.media_db = null
params.media = null
params.method = "carveme"
params.threads = 12


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
      --method [str]                The algorithm to use. Either `carveme` or `gapseq`. `gapseq`
                                    requires docker or singularity.
    Growth Media:
      --media_db                    A file containing growth media specification.
                                    `*.tsv` for CARVEME and `*.csv` for gapseq.
      --media                       Comma-separated list of media names to use. Only used for CARVEME.
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
  from os import path
  from carveme import config, project_dir
  from carveme.cli.carve import first_run_check

  diamond_db = project_dir + config.get('generated', 'diamond_db')[:-5] + ".dmnd"

  if __name__ == "__main__":
    if path.exists(diamond_db):
      subprocess.check_output(["rm", diamond_db])
    first_run_check()
    with open("db_stats.txt", "w") as out:
      res = subprocess.check_output(["diamond", "dbinfo", "-d", diamond_db])
      out.write(res.decode())
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

process checkm {
  cpus params.threads
  publishDir "${params.data_dir}", mode: "copy", overwrite: true

  input:
  path(proteins)

  output:
  path("checkm_summary.tsv")

  """
  checkm lineage_wf --genes -t ${task.cpus} -x faa . checkm
  checkm qa checkm/lineage.ms checkm --tab_table -f checkm_summary.tsv
  """
}

process build_carveme {
  cpus 2
  publishDir "${params.data_dir}/carveme_models", mode: "copy", overwrite: true

  input:
  tuple val(id), path(genes_dna), path(genes_aa), path(db_info)

  output:
  tuple val("${id}"), path("${id}.xml.gz"), path("${id}.log")

  script:
  if (params.media_db && params.media)
    """
    CPX_PARAM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} \
    carve ${genes_aa} -o ${id}.xml.gz --mediadb ${params.media_db} \
          --gapfill ${params.media} --diamond-args "-p ${task.cpus}" \
          --fbc2 -v > ${id}.log
    """
  else if (params.media)
    """
    CPX_PARAM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} \
    carve ${genes_aa} -o ${id}.xml.gz --gapfill ${params.media} \
          --diamond-args "-p ${task.cpus}" \
          --fbc2 -v > ${id}.log
    """
  else
    """
    CPX_PARAM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} \
    carve ${genes_aa} -o ${id}.xml.gz --diamond-args "-p ${task.cpus}" \
      --fbc2 -v > ${id}.log
    """
}

process build_gapseq {
  cpus 1
  publishDir "${params.data_dir}/gapseq_models", mode: "copy", overwrite: true

  input:
  tuple val(id), path(assembly)

  output:
  tuple val("${id}"), path("${id}.xml.gz"), path("${id}.log")

  script:
  if (params.media_db)
    """
    cp ${baseDir}/${params.media_db} medium.csv
    gapseq doall ${assembly} medium.csv > ${id}.log
    gzip ${id}.xml
    """
  else
    """
    gapseq doall ${assembly} /opt/gapseq/dat/media/gut.csv > ${id}.log
    gzip ${id}.xml
    """
}

process check_model {
  cpus 1
  publishDir "${params.data_dir}/model_qualities", mode: "copy", overwrite: true

  input:
  tuple val(id), path(model), path(log)

  output:
  tuple val("${id}"), path("${id}.html.gz")

  script:
  if (params.method == "gapseq")
    """
    CPX_PARAM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} \
    memote report snapshot ${model} --filename ${id}.html
    gzip ${id}.html
    """
  else
    """
    CPX_PARAM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} \
    memote report snapshot ${model} --solver cplex --filename ${id}.html
    gzip ${id}.html
    """
}

process carveme_fba {
  cpus 1

  input:
  tuple val(id), path(model), path(log)

  output:
  tuple val(id), path("${id}_exchanges.csv"), path("${id}_growth_rate.csv")

  """
  #!/usr/bin/env python3

  import cobra
  import pandas as pd
  from carveme import project_dir
  from os import path

  model = cobra.io.read_sbml_model("${model}")

  exids = [r.id for r in model.exchanges]
  if "${params.media}" != "null":
    if "${params.media_db}" == "null":
      media_df = pd.read_csv(
        path.join(project_dir, "${model}", "input", "media_db.tsv"), sep="\\t")
    else:
      media_df = pd.read_csv("${params.media_db}", sep="\\t")
    mname = "${params.media}".split(",")[0]
    media_df = media_df[media_df.medium == mname]
    if "flux" not in media_df.columns:
      media_df["flux"] = 0.1
    if "reaction" not in media_df.columns:
      media_df["reaction"] = "EX_" + media_df["compound"] + "_e"
    media_df.index = media_df.reaction
    model.medium = media_df.flux[media_df.index.isin(exids)]

  rate = pd.DataFrame({"id": "${id}", "growth_rate": model.slim_optimize()}, index=[0])
  sol = cobra.flux_analysis.pfba(model)
  ex_fluxes = sol.fluxes[
    sol.fluxes.index.isin(exids)
    & (sol.fluxes.abs() > model.tolerance)
  ]
  met_names = ex_fluxes.index.to_series().apply(
    lambda idx: model.metabolites.get_by_id(idx.replace("EX_", "")).name)
  exchanges = pd.DataFrame({
    "assembly": "${id}",
    "reaction": ex_fluxes.index,
    "flux": ex_fluxes,
    "description": met_names
  })
  rate.to_csv("${id}_growth_rate.csv", index=False)
  exchanges.to_csv("${id}_exchanges.csv", index=False)
  """
}

process summarize_fba {
  cpus params.threads
  publishDir "${params.data_dir}", mode: "copy", overwrite: true

  input:
  path(exchanges)
  path(rates)

  output:
  tuple path("exchanges.csv"), path("growth_rates.csv")

  """
  #!/usr/bin/env python3

  import pandas as pd
  import glob

  res = map(pd.read_csv, glob.glob("*_exchanges.csv"))
  exchanges = pd.concat(res)
  exchanges.to_csv("exchanges.csv", index=False)

  res = map(pd.read_csv, glob.glob("*_growth_rate.csv"))
  growth_rates = pd.concat(list(res))
  growth_rates.to_csv("growth_rates.csv", index=False)
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

  def models = null
  if (params.method == "carveme") {
    init_db()
    find_genes(genomes)
    build_carveme(find_genes.out.combine(init_db.out))
    models = build_carveme.out
    models | carveme_fba
    exchanges = carveme_fba.out.map{it -> it[1]}.collect()
    rates = carveme_fba.out.map{it -> it[2]}.collect()
    summarize_fba(exchanges, rates)
    checkm(find_genes.out.map{prots -> prots[2]}.collect())
  } else if (params.method == "gapseq") {
    build_gapseq(genomes)
    models = build_gapseq.out
  } else {
    error "Method must be either `carveme` or `gapseq`."
  }

  check_model(models)
}
