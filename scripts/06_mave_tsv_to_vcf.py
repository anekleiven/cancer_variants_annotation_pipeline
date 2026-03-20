
#====================================================================
# TSV to VCF formatting 
# (Preparation for MAVE annotation)
# ====================================================================

""" 
Script: 06_mave_tsv_to_vcf.py 
Author: Ane Kleiven 

Major outputs: 
  Prepare variant data for assembly liftover and MAVE annotation
  Converting .tsv files to .vcf format 

"""

# -------------------------------------------
# Import libraries 
# -------------------------------------------

import argparse
from pathlib import Path
import pandas as pd
import numpy as np 

# -------------------------------------------
# Argparse function for user input file paths
# -------------------------------------------

def getargs(): 
  parser = argparse.ArgumentParser(
    description="Prepare variant data for MAVE annotation. Formatting from .tsv to .vcf"
  ) 

  parser.add_argument(
    "--somatic_variants", 
    type=Path, 
    required=False,
    default="output/variants_with_germline_proximity.tsv",
    help="Path to the variant input file (e.g. output/variants_with_germline_proximity.tsv)"
  )

  parser.add_argument(
     "--neutral_clinvar",
     type=Path, 
     required=False,
     default= "output/neutral_clinvar_filtered.tsv",
     help="Path to the neutral variant file from ClinVar"
  )

  parser.add_argument(
    "--somatic_output",
    type=Path,
    required=False,
    default="output/somatic_variants_GRCh37.vcf",
    help="Path to the annotated output variant file (e.g. output/somatic_variants_GRCh37.vcf)"
  )

  parser.add_argument(
     "--clinvar_output",
     type=Path, 
     required=False,
     default= "output/neutral_clinvar_GRCh37.vcf",
     help="Path to the Clinvar vcf output file."
  )

  return parser.parse_args() 

# -------------------------------------------
# TSV to VCF function 
# -------------------------------------------

def tsv_to_vcf(df, mapping, output_path, source_name="Ane_Kleiven_Pipeline"):
  """
  Converts .tsv variant files to .vcf format
  Mapping: dictionary which maps vcf entries to df columns
  Example: {'CHROM':'Chromosome', 'POS':'Position}
  """

  print("-"*50)
  print("TSV TO VCF CONVERTER🧬")
  print("-"*50 + "\n") 

  initial_count = len(df)
  print(f"Converting {initial_count} variants to VCF..")

  # Mapping ID
  if 'ID' in mapping and mapping['ID'] in df.columns:
      ids = df[mapping['ID']].astype(str)
  else:
      ids = "."

  # Create vcf df 
  vcf_df = pd.DataFrame({
      "#CHROM": "chr" + df[mapping['CHROM']].astype(str).str.replace('chr', ''),
      "POS": df[mapping['POS']].astype(int),
      "ID": ids,
      "REF": df[mapping['REF']],
      "ALT": df[mapping['ALT']],
      "QUAL": ".",
      "FILTER": "PASS",
      "INFO": "GENE=" + df[mapping['GENE']].astype(str)
  })
  
  # Add EXTRA_INFO if accessible (ClinicalSignificance, ReviewStatus etc.)
  if 'EXTRA_INFO' in mapping:
      for key, col in mapping['EXTRA_INFO'].items():
          vcf_df["INFO"] += f";{key}=" + df[col].astype(str).str.replace(" ", "_")

  # Clean data 
  vcf_df = vcf_df.dropna(subset=["REF", "ALT"])
  vcf_df = vcf_df[(vcf_df["REF"] != "-") & (vcf_df["ALT"] != "-")]
  vcf_df = vcf_df[vcf_df["REF"] != vcf_df["ALT"]]
  vcf_df = vcf_df.sort_values(by=["#CHROM", "POS"], ascending=True)
  vcf_df = vcf_df[(vcf_df["REF"].str.len() == 1) & (vcf_df["ALT"].str.len() == 1)]    # For simplicity: keep only SNVs 
  
  final_count = len(vcf_df)
  variants_removed = initial_count - final_count

  print(f"Total variants in: {initial_count}")
  print(f"Variants removed: {variants_removed}")
  print(f"Final variant count: {final_count}\n")

  # Write file with header
  vcf_header = [
      "##fileformat=VCFv4.2",
      f"##source={source_name}",
      "##INFO=<ID=GENE,Number=1,Type=String,Description='Gene Symbol'>",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  ]

  with open(output_path, "w") as f:
      for line in vcf_header:
          f.write(line + "\n")
      vcf_df.to_csv(f, sep="\t", index=False, header=False)
  
  print(f"VCF saved to {output_path}.\n")


# -------------------------------------------
# Main function 
# -------------------------------------------

def main():
  
  args = getargs() 

  # Convert somatic variant file to VCF
  somatic_df = pd.read_csv(args.somatic_variants, sep="\t", low_memory=False)
  somatic_mapping = {
    "CHROM": "Chromosome",
    "POS": "Start_Position",
    "ID": ".",
    "REF": "Reference_Allele",
    "ALT": "Tumor_Seq_Allele2",
    "GENE": "Hugo_Symbol"}
  
  tsv_to_vcf(somatic_df, somatic_mapping, args.somatic_output)

  # Convert neutral ClinVar file to VCF 
  clinvar_df = pd.read_csv(args.neutral_clinvar, sep="\t", low_memory=False) 
  clinvar_mapping = {
     "CHROM": "Chromosome",
     "POS": "PositionVCF",
     "ID": "VariationID",
     "REF": "ReferenceAlleleVCF",
     "ALT": "AlternateAlleleVCF",
     "GENE":"GeneSymbol", 
     "EXTRA_INFO": {"CLINSIG": "ClinicalSignificance", "REVSTAT": "ReviewStatus"}
    }
  
  tsv_to_vcf(clinvar_df, clinvar_mapping, args.clinvar_output)

 
if __name__ == "__main__":
   main() 