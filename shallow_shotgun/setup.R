#! /usr/bin/env Rscript

install.packages("BiocManager")
remotes::install_github("gibbons-lab/mbtools", repos = BiocManager::repositories())
