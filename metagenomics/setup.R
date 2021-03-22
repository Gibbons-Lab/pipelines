#! /usr/bin/env Rscript

install.packages(c("BiocManager", "remotes"))
remotes::install_github("gibbons-lab/mbtools", repos = BiocManager::repositories())
