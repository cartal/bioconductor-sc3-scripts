#!/usr/bin/env Rscript 
# Script to creare SCE object for clustering with SC3.
# Carlos Talavera-López - 2018

# Load required packages
packages <- c("tidyverse", "SingleCellExperiment", "scater", "optparse", "SC3","workflowscriptscommon")
suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))

# parse options

option_list = list(
  make_option(
    c("-o", "--SCE"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A QC-ed, nornmalised SingleCellExperiment object from scater."
  ),
  make_option(
    c("-g", "--gene_filter"),
    action = 'store',
    default = NA,
    type = 'character',
    help = "A boolean variable which defines whether to perform gene filtering before SC3 clustering."
  ),
  make_option(
    c("-pdmin", "--pct_dropout_min"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "If gene_filter = TRUE, then genes with percent of dropouts smaller than pct_dropout_min are filtered out before clustering."
  ),
  make_option(
    c("-pdmax", "--pct_dropout_max"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "If gene_filter = TRUE, then genes with percent of dropouts larger than pct_dropout_max are filtered out before clustering."
  ),
  make_option(
    c("-drmin", "--d_region_min"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "Defines the minimum number of eigenvectors used for kmeans clustering as a fraction of the total number of cells. Default is 0.04. See SC3 paper for more details."
  ),
  make_option(
    c("-drmax", "--d_region_max"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "Defines the maximum number of eigenvectors used for kmeans clustering as a fraction of the total number of cells. Default is 0.07. See SC3 paper for more details."
  ),
  make_option(
    c("-svm_num", "--svm_num_cells"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "Number of randomly selected training cells to be used for SVM prediction. Default is NULL."
  ),
  make_option(
    c("-svm_i", "--svm_train_indeces"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "A numeric vector defining indeces of training cells that should be used for SVM training. The default is NULL."
  ),
  make_option(
    c("-svm_max", "--svm_max_cells"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "Define the maximum number of cells below which SVM are not run."
  ),
  make_option(
    c("-t", "--num_threads"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "Number of threads/cores to be used in the user's machine."
  ),
  make_option(
    c("-seed", "--rand_seed"),
    action = "store",
    default = 8,
    type = 'character',
    help = "RDS output object"
  ),
  make_option(
    c("-o", "--output"),
    action = "store",
    default = NA,
    type = 'character',
    help = "RDS output object"
  )
)


opt <- wsc_parse_args(option_list)

# Check parameter values

if ( ! dir.exists(opt$SCE)){
  stop((paste('The SCE object', opt$SCE, 'does not exist')))
}

# Create SCE object from data

sc3_obj <- sc3_prepare(opt$SCE, gene_filter = opt$gene_filter, opt$pct_dropout_min = 10, opt$pct_dropout_max = 90, opt$d_region_min = 0.04,
  opt$d_region_max = 0.07, opt$svm_num_cells = NULL, opt$svm_train_indeces = NULL, opt$svm_max_cells = 5000, opt$num_threads = 8, kmeans_nstart = NULL,
  kmeans_iter_max = 1e+09, rand_seed = 1712)

# Output to a serialized R object

saveRDS(sc3_obj, file = opt$output)
