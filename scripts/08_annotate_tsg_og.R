
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
  library(tidyverse)
  library(stringr)
  library(data.table)
})

# -----------------------------
# Parse command-line arguments
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  message("\nUsage: Rscript scripts/08_annotate_tsg_og.R <input_file> <output_file> [cache_dir]")
  message("Example: Rscript scripts/08_annotate_tsg_og.R output/variants_tsg_og.tsv data_cache\n")
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

variants <- data.table::fread(input_file, 
                  sep = "\t", 
                  colClasses = "character", 
                  na.strings = c("", "NA"))

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
  mutate(Entrez_Gene_Id = as.numeric(Entrez_Gene_Id))


variants <- variants |>
  dplyr::left_join(gene_annotations_filtered, by = c("Entrez_Gene_Id" = "entrezgene"))

# view annotated data 
variants |>
  dplyr::select(Entrez_Gene_Id, Hugo_Symbol, ncg_tsg, ncg_oncogene) |>
  head(20)


# -----------------------------
# Create col 'is_null_var_tsg' 
# -----------------------------

# check values in "Variant_Type"
table(variants$Consequence)

# define consequence types for null variants 
null_pattern <- paste(
  "stop_gained",                   # Nonsense
  "frameshift_variant",            # Frameshift
  "splice_acceptor_variant",       # Canonical splice ±1,2
  "splice_donor_variant",          # Canonical splice ±1,2
  "start_lost",                    # Initiation codon
  sep = "|"
)

# set NA values in ncg_tsg and ncg_oncogene to FALSE
variants <- variants |>
  mutate(
    ncg_tsg = as.logical(replace_na(ncg_tsg, FALSE)),
    ncg_oncogene = as.logical(replace_na(ncg_oncogene, FALSE))
  )

# create null variant  and null variant in tsg column 
variants <- variants |>
  mutate(
    is_null_variant = str_detect(Consequence, null_pattern),
    is_null_var_tsg = is_null_variant & (ncg_tsg == TRUE)
  ) |>
  mutate(
    is_null_variant = replace_na(is_null_variant, FALSE),
    is_null_var_tsg = replace_na(is_null_var_tsg, FALSE)
  )|>
  mutate(
    across(everything(), as.character)
  )

# view results
message("Number of null variants in dataset:")
table(variants$is_null_variant)
message("\nNumber of null variants in tumor suppressor genes:")
table(variants$is_null_var_tsg)

# Check oncogenicity 
message("\nCross-tab: All Null Variants vs Oncogenicity")
table(Oncogenic = variants$ONCOGENIC, Null_Variant = variants$is_null_variant)

# Check oncogenicity 
message("\nCross-tab: Null Variants in TSG vs Oncogenicity")
table(Oncogenic = variants$ONCOGENIC, Null_in_TSG = variants$is_null_var_tsg)

# -----------------------------
# Save annotate dataset to .tsv 
# -----------------------------

variants_final <- variants %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~replace_na(., "")))

message("\nSaving annotated variants to file...")
fwrite(variants_final, 
       output_file, 
       sep = "\t",      
       na = "",        
       quote = FALSE,
       scipen = 999)

message("\nDone! Variants with TSG and OG annotations written to:")
message(output_file)
message("-----------------------------------------------------------\n")