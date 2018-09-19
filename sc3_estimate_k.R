#!/usr/bin/env Rscript 
# Uses Tracy-Widom theory on random matrices to estimate the optimal number of clusters k. It creates and populates the k_estimation item of the sc3 slot of the metadata(object)..
# Carlos Talavera-López - 2018

# Load required packages
packages <- c("tidyverse", "SingleCellExperiment", "optparse", "SC3","workflowscriptscommon")
suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))

# parse options

option_list = list(
  make_option(
    c("-o", "--SCE"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A Single Cell Experiment object previously prepared with SC3."
  ),
  make_option(
    c("-out", "--output"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Output k"
  )
)


opt <- wsc_parse_args(option_list)

# Check if all parameters have been provided by the user

if ( ! dir.exists(opt$SCE)){
  stop((paste('The SCE object', opt$SCE, 'does not exist')))
}

# Select K 

sc3_k <- sc3_estimate_k(opt$SCE)

# Output to a serialized R object

saveRDS(sc3_k, file = opt$out)
