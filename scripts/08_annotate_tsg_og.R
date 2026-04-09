
#====================================================================
# Tumor suppressor- and oncogene annotation 
#====================================================================


# Script: 08_annotate_tsg_og.R
# Author: Ane Kleiven

# Description:
  # Annotates a variant dataset with tumor suppressor gene- and oncogene annotations 
  # using the geneOncoX package (Sigve Nakken)

# Usage:
  # Rscript 08_annotate_tsg_og.R <input_file> <output_file> [cache_dir]

# -----------------------------
# Load required libraries
# -----------------------------

suppressPackageStartupMessages({
  library(tidyverse)})

# -----------------------------
# Parse command-line arguments
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  message("\nUsage: Rscript scripts/08_annotate_tsg_og.R <input_file> <output_file> [cache_dir]")
  message("Example: Rscript scripts/08_annotate_tsg_og.R output/variants_with_maves.tsv data/variants_tsg_og.tsv data_cache\n")
  stop("Missing required arguments. Please provide at least input and output file paths.")
}

input_file  <- args[1]
output_file <- args[2]
cache_dir   <- ifelse(length(args) >= 3, args[3], "cancer_variants_annotation_pipeline/data_cache")

# -----------------------------
# Display settings
# -----------------------------
message("-----------------------------------------------------------")
message(" UniProt Mapping Script ")
message("-----------------------------------------------------------")
message("Input file:  ", input_file)
message("Output file: ", output_file)
message("Cache dir:   ", cache_dir)
message("-----------------------------------------------------------\n")

# -----------------------------
# Ensure cache directory exists
# -----------------------------
if (!dir.exists(cache_dir)) {
  message("\nCreating cache directory: ", cache_dir)
  dir.create(cache_dir, recursive = TRUE)
}

# -----------------------------
# Install geneOncoX if missing
# -----------------------------
if (!requireNamespace("geneOncoX", quietly = TRUE)) {
  message("\nInstalling 'geneOncoX' from GitHub (first time only)...")
  remotes::install_github("sigven/geneOncoX", quiet = TRUE)
}
# -----------------------------
# Download Get basic
# cancer-relevant gene annotations
# -----------------------------

message("Loading get_basic dataset from geneOncoX...\n")
gene_annotations <- geneOncoX::get_basic(cache_dir = cache_dir, force_download = FALSE)

# -----------------------------
# Load variant dataset
# -----------------------------
message("\nReading variant data..")
variants <- readr::read_tsv(input_file, show_col_types = FALSE)
message("\nLoaded ", nrow(variants), " variants.\n")


# -----------------------------
# Annotate dataset with
# TSG and OG annotations 
# -----------------------------

# filter gene annotations to entrez gene id, oncogene and tumor suppressor gene annotations 
gene_annotations_filtered <- gene_annotations$records |>
  dplyr::select(entrezgene, ncg_tsg, ncg_oncogene)

head(gene_annotations_filtered) 

# join tables 
variants <- variants |>
  dplyr::left_join(gene_annotations_filtered, by = c("Entrez_Gene_Id" = "entrezgene"))

# view annotated data 
variants |>
  dplyr::select(Entrez_Gene_Id, Hugo_Symbol, ncg_tsg, ncg_oncogene) |>
  head(20)


# -----------------------------
# Save annotate dataset to .tsv 
# -----------------------------

message("\nSaving annotated variants to file...")
readr::write_tsv(variants, output_file)

message("\nDone! Variants with TSG and OG annotations written to:")
message(output_file)
message("-----------------------------------------------------------\n")