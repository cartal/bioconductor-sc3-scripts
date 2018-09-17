#!/usr/bin/env Rscript 
# Script to creare SCE object for clustering with SC3.
# Carlos Talavera-López - 2018

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Load common functions

suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options

option_list = list(
  make_option(
    c("-c", "--counts"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A matrix of expression data with rows corresponding to genes and columns corresponding to cells."
  ),
  make_option(
    c("-a", "--annotation"),
    action = "store",
    default = NA,
    type = 'character',
    help = "An annotation matrix where rows are the cell names (column names in the expression matrix)"
  ),
  make_option(
    c("-o", "--output"),
    action = "store",
    default = NA,
    type = 'character',
    help = "RDS output object"
  ),
)

opt <- wsc_parse_args(option_list, mandatory = c('cts', 'ann', 'out'))

# Check parameter values

if ( ! dir.exists(opt$cts)){
  stop((paste('The matrix', opt$cts, 'does not exist')))
}

# Load SC3 and dependencies

packages <- c("tidyverse", "SC3", "SingleCellExperiment", "scater")
suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))

# Create SCE object from data

sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(opt$cts),
        logcounts = log2(as.matrix(opt$cts) + 1)
    ), 
    colData = opt$ann
)

# Add gene/transcript symbol

rowData(sce)$feature_symbol <- rownames(sce)

# Remove duplicated entries

sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

# Define Spike-in's 

isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)

# Output to a serialized R object

saveRDS(sce, file = opt$out)
