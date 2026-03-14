# -----------------------------------------------------------
#====================================================================
# Uniprot mapping script
#====================================================================

# Script: 02_add_uniprot_mapping.R
# Author: Ane Kleiven

# Description:
#   Annotates a variant dataset with UniProt accession numbers
#   using the geneOncoX reference (Sigve Nakken)
#
# Usage:
#   Rscript 02_add_uniprot_mapping.R <input_file> <output_file> [cache_dir]
#
# Example:
#   Rscript scripts/02_add_uniprot_mapping.R \
#     output/annotated_with_hotspots.maf \
#     data/variants_with_uniprot.tsv \
#     data_cache
# -----------------------------------------------------------

# -----------------------------
# Load required libraries
# -----------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(remotes)
})

# -----------------------------
# Parse command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  message("\nUsage: Rscript scripts/02_add_uniprot_mapping.R <input_file> <output_file> [cache_dir]")
  message("Example: Rscript scripts/02_add_uniprot_mapping.R output/annotated_with_hotspots.maf data/variants_with_uniprot.tsv data_cache\n")
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
# Load gene reference data
# -----------------------------
message("Loading gene reference data from geneOncoX...\n")
gene_xref <- geneOncoX::get_gencode(cache_dir = cache_dir)


# -----------------------------
# Prepare gene–UniProt mapping
# -----------------------------
message("\nPreparing gene symbol to UniProt mapping table...")
gene_to_uniprot <- gene_xref$records$grch37 %>%
  dplyr::select(symbol, uniprot_acc) %>%
  dplyr::filter(!is.na(uniprot_acc)) %>%
  dplyr::distinct()

message("Loaded ", nrow(gene_to_uniprot), " unique gene to UniProt mappings.\n")

# -----------------------------
# Load variant dataset
# -----------------------------
message("\nReading variant data..")
variants <- readr::read_tsv(input_file, show_col_types = FALSE)
message("\nLoaded ", nrow(variants), " variants.\n")

# -----------------------------
# Map gene symbols to UniProt accessions
# -----------------------------
message("\nMapping variants to UniProt accessions..\n")
variants_with_uniprot <- variants %>%
  dplyr::left_join(gene_to_uniprot, by = c("Hugo_Symbol" = "symbol"))

# separate mapped and unmapped
mapped <- variants_with_uniprot %>% dplyr::filter(!is.na(uniprot_acc))
unmapped  <- variants_with_uniprot %>% dplyr::filter(is.na(uniprot_acc)) %>% dplyr::select(-uniprot_acc)

unmapped_genes <- variants_with_uniprot %>%
  dplyr::filter(is.na(uniprot_acc)) %>%
  count(Hugo_Symbol, sort = TRUE)

head(unmapped_genes, 10)

# count mapped and unmapped 
unmapped_count <- sum(is.na(variants_with_uniprot$uniprot_acc))
message("\nFinal Mapping Results:")
message("Total variants:    ", nrow(variants_with_uniprot))
message("Successfully mapped: ", nrow(variants_with_uniprot) - unmapped_count)
message("Unmapped:      ", unmapped_count)

# -----------------------------
# Save results
# -----------------------------

message("\nSaving annotated variants to file...")
readr::write_tsv(variants_with_uniprot, output_file)

message("\nDone! Variants with UniProt accessions written to:")
message(output_file)
message("-----------------------------------------------------------\n")