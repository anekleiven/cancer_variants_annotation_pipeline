""" 
====================================================================
Germline Proximity Annotation Script
====================================================================

Script: 06_annotate_germline_proximity.py 
Author: Ane Kleiven 

Major outputs: 
  1. Loading variant data files 
  2. Preprocessing of germline variant file (filter variants, extract AA position) 
     Can be found here: (https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/)
  3. Create dictionary of germline variants per gene, with AA positions as values
  4. Calculate the shortest germline proximity for each somatic variant 
  5. Apply germline proximity to somatic variant table 
  7. Save new output file 

"""

# import libraries 

import argparse
from pathlib import Path
import pandas as pd
import numpy as np 

# create argparse function for user input file paths 

def getargs(): 
  parser = argparse.ArgumentParser(
    description="Annotate somatic cancer variants with germline variant proximity"
  ) 

  parser.add_argument(
    "--input", 
    type=Path, 
    default="annotation_pipeline/output/variants_with_func_sites.tsv",
    help="Path to the variant input file (e.g. annotation_pipeline/output/variants_with_func_sites.tsv)"
  )

  parser.add_argument(
    "--germline_variant_file",
    type=Path,
    default="annotation_pipeline/data/variant_summary.txt.gz",
    help="Path to the germline variant file (e.g. annotation_pipeline/data/variant_summary.txt)"
  )

  parser.add_argument(
    "--output",
    type=Path,
    default="annotation_pipeline/output/variants_with_germline_proximity.tsv",
    help="Path to the annotated output variant file (e.g. annotation_pipeline/output/variants_with_germline_proximity.tsv)"
  )

  return parser.parse_args() 



# create main function 
def main(): 

  args = getargs() 

  # load files needed for annotation 
  print("Loading files needed for annotation of germline proximity..\n")

  somatic_variants = pd.read_csv(args.input, sep="\t", low_memory=False) 
  germline_variants=pd.read_csv(args.germline_variant_file, sep="\t", low_memory=False)

  print(f"Loaded {len(somatic_variants):,} somatic variants from {args.input.name}") 
  print(f"Loaded {len(germline_variants):,} variants from {args.germline_variant_file.name}\n")

  # Extract variants that are:
  #   From assembly GRCh38
  #   Pathogenic and Likely Pathogenic variants 
  #   Origin == germline

  print("Filtering variants to only include germline Pathogenic/LikelyPathogenic variants from assembly GRCh38..\n")

  germline_variants_filtered = germline_variants[
      (germline_variants["Assembly"] == "GRCh38") &
      (germline_variants["Origin"] == "germline") &
      (germline_variants["ClinicalSignificance"].isin(["Pathogenic", "LikelyPathogenic"]))
  ]

  print(f"Number of germline variants after filtering: {len(germline_variants_filtered):,}\n") 

  # Extract amino acid position (p.Gly56Arg) == position 56 
  # Exluding frameshift and complex variants

  print("Extracting HGVSp from 'Name' column..\n")

  germline_variants_filtered["HGVSp"] = (
      germline_variants_filtered["Name"]
      .str.split(" ", n=1)
      .str[1]
      .str.strip("()")
  )

  print("Extracting amino acid position from HGVSp using regex..\n")

  germline_variants_filtered["AA_Position_Germline"] = (
    germline_variants_filtered["HGVSp"]
    .str.extract(r"p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}")
    .astype(float)
  )

  print("Preview of the 'AA_Position_Germline' column:")
  print(germline_variants_filtered["AA_Position_Germline"].head(), "\n")


  # Clean germline variant data before mapping

  print("Keep only selected columns for further analysis..\n")

  columns_mapping = ["GeneSymbol",
                     "GeneID",
                     "HGVSp",
                     "AA_Position_Germline",
                     "ClinicalSignificance"
                     ]
  
  germline_variants_cleaned = germline_variants_filtered[columns_mapping]

  print("Remove rows with missing data..\n")

  germline_variants_cleaned = germline_variants_cleaned.dropna()

  print("Rename GeneID column to match somatic variant file..\n")

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
  print(dict(list(germline_per_gene.items())[:3]), "\n")      # first three key-value pairs

  
  # define a function to calculate germline proximity per somatic variant (row by row)

  def compute_germline_proximity(row):
      gene_id = row["Entrez_Gene_Id"]
      somatic_pos = row["Protein_position"]

      if gene_id not in germline_per_gene or pd.isna(somatic_pos):
          return np.nan

      germline_positions = germline_per_gene[gene_id]
      return min(abs(somatic_pos - gp) for gp in germline_positions)


  # call the function 

  print("Computing germline proximity and applying distances to somatic variants..\n")

  somatic_variants["Germline_Proximity"] = somatic_variants.apply(
    compute_germline_proximity,
    axis=1
)

  # look at the data after merging: 
  selected_columns = ["Hugo_Symbol", "Entrez_Gene_Id", "HGVSp", "Protein_position", "Germline_Proximity"]

  print("Somatic variant data after merging:")
  print(somatic_variants[selected_columns].head(), "\n")


  # save final annotated file: 
  somatic_variants.to_csv(args.output, sep="\t", index=False)
  print(f"Annotated output file saved as: {args.output}\n")


if __name__ == "__main__":
	main() 
