
# ====================================================================
# Germline Proximity Annotation 
# ====================================================================

""" 
Script: 05_annotate_germline_proximity.py 
Author: Ane Kleiven 

This script requires the ClinVar variant file. 
The ClinVar file can be accessed here: 
     (https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/)

Major outputs: 
  1. Loading variant data files
  2. Preprocessing of ClinVar variant file (filter variants, extract AA position) 
  3. Create dictionary of germline variants per gene, with AA positions as values
  4. Calculate the shortest germline distance for each somatic variant 
  5. Apply germline distances to somatic variant table 
  7. Save new output file (.tsv) 

"""

#-----------------------------------------------------
# Import libraries 
#-----------------------------------------------------

import argparse
from pathlib import Path
import pandas as pd
import numpy as np 

#-----------------------------------------------------
# Create argparse function for user input file paths 
#-----------------------------------------------------

def getargs(): 
  parser = argparse.ArgumentParser(
    description="Annotate somatic cancer variants with distance to known pathogenic germline variants."
  ) 

  parser.add_argument(
    "--input", 
    type=Path, 
    required=False,
    default="output/variants_with_func_sites.tsv",
    help="Path to the variant input file (e.g. output/variants_with_func_sites.tsv)"
  )

  parser.add_argument(
    "--germline_variant_file",
    type=Path,
    required=False, 
    default="data/variant_summary.txt.gz",
    help="Path to the germline variant file (e.g. data/variant_summary.txt)"
  )

  parser.add_argument(
    "--output",
    type=Path,
    required=False,
    default="output/variants_with_germline_proximity.tsv",
    help="Path to the annotated output variant file (e.g. output/variants_with_germline_proximity.tsv)"
  )

  return parser.parse_args() 


#-----------------------------------------------------
# Main() function 
#-----------------------------------------------------

def main(): 

  args = getargs() 

  # Load files needed for annotation 
  print("Loading files needed for annotation of germline proximity..\n")

  print("Reading somatic variant file..")
  somatic_variants = pd.read_csv(args.input, sep="\t", low_memory=False) 
  print(f"Loaded {len(somatic_variants):,} somatic variants from {args.input.name}\n")

  print("Reading germline variant file (filtered)..")
  # Choose germline columns to read 
  germline_cols = [
    "Assembly",
    "Origin",
    "Name",
    "GeneID",
    "GeneSymbol", 
    "ClinicalSignificance"
  ]

  # Only include gene ID's present in the somatic variant file
  somatic_gene_ids = set(somatic_variants["Entrez_Gene_Id"].dropna().unique())

  # Read file in chunks to save memory 
  chunks =  []

  for chunk in pd.read_csv(
    args.germline_variant_file, 
    sep="\t", 
    usecols=germline_cols, 
    chunksize=100000,
    low_memory=False
  ):
    chunks.append(chunk[chunk["GeneID"].isin(somatic_gene_ids)])
  
  germline_variants = pd.concat(chunks, ignore_index=True)

  print(f"Loaded {len(germline_variants):,} variants from {args.germline_variant_file.name}\n")

  # Extract variants that are:
  #   From assembly GRCh37 (Assembly)
  #   Pathogenic and Likely Pathogenic variants (ClinicalSignificance)
  #   Germling (Origin) 

  print("Filtering variants to only include germline Pathogenic/Likely Pathogenic variants from assembly GRCh37..\n")

  germline_variants_filtered = germline_variants[
      (germline_variants["Assembly"] == "GRCh37") &
      (germline_variants["Origin"] == "germline") &
      (germline_variants["ClinicalSignificance"].isin(["Pathogenic", "LikelyPathogenic"]))
  ]

  print(f"Number of germline variants after filtering: {len(germline_variants_filtered):,}\n") 


  # Extract amino acid position ((p.Gly56Arg) == position 56)
  # Exluding frameshift and complex variants

  print("Extracting HGVSp from 'Name' column..\n")
  germline_variants_filtered["HGVSp"] = (
      germline_variants_filtered["Name"]
      .str.split(" ", n=1)
      .str[1]
      .str.strip("()")
  ).copy()

  print("Extracting amino acid position from HGVSp using regex..\n")
  germline_variants_filtered["AA_Position_Germline"] = (
    germline_variants_filtered["HGVSp"]
    .str.extract(r"p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}")
    .astype(float)
  )

  print("Preview of the 'AA_Position_Germline' column:")
  print(germline_variants_filtered["AA_Position_Germline"].head(), "\n")


  # Clean germline variant data before mapping
  print("Cleaning germline variant data before mapping...")

  print("Removing rows with missing data..\n")
  germline_variants_cleaned = germline_variants_filtered.dropna()

  print("Renaming columns...\n")
  germline_variants_cleaned = germline_variants_cleaned.rename(
      columns={"GeneID": "Entrez_Gene_Id",
               "HGVSp":"HGVSp_germline"}
  )

  print("Germline variant file after cleaning:")
  print(germline_variants_cleaned.head(), "\n")


  # Create dictionary with genes as keys and AA positions as values
  print("Grouping germline variants after genes with AA positions into lists..\n")
  germline_per_gene = (
      germline_variants_cleaned
      .groupby("Entrez_Gene_Id")["AA_Position_Germline"]
      .apply(list)
      .to_dict()
  )

  print("Preview of the grouped germline dictionary:")
  print(dict(list(germline_per_gene.items())[:5]), "\n")    

  #-----------------------------------------------------
  # Calculate germline proximity per somatic variant 
  # (row by row)
  #-----------------------------------------------------

  def compute_germline_proximity(row):
      gene_id = row["Entrez_Gene_Id"]
      somatic_pos = row["Protein_position"]

      if gene_id not in germline_per_gene or pd.isna(somatic_pos):
          return np.nan

      germline_positions = germline_per_gene[gene_id]
      return min(abs(somatic_pos - gp) for gp in germline_positions)


  # Call germline proximity function
  print("Computing germline proximity and applying distances to somatic variants..\n")
  somatic_variants["Germline_Proximity"] = somatic_variants.apply(
    compute_germline_proximity,
    axis=1
)

  print("Somatic variant data after merging:")
  selected_columns = ["Hugo_Symbol", "Entrez_Gene_Id", "HGVSp", "Protein_position", "Germline_Proximity"]
  print(somatic_variants[selected_columns].head(), "\n")


  # Save final annotated file: 
  somatic_variants.to_csv(args.output, sep="\t", index=False)
  print(f"Annotated output file saved as: {args.output}\n")


  print("Script ended!🎉\n")


if __name__ == "__main__":
	main() 
