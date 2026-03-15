""" 
====================================================================
Somatic Variant file (.tsv) to VCF format 
Preparation for MAVE annotation 
====================================================================

Script: 06_mave_tsv_to_vcf.py 
Author: Ane Kleiven 

Major outputs: 
  Prepare variant data for MAVE annotation using VEP
  Converting .tsv file to .vcf format 

"""

# import libraries 
import argparse
from pathlib import Path
import pandas as pd
import numpy as np 


# create argparse function for user input file paths 
def getargs(): 
  parser = argparse.ArgumentParser(
    description="Prepare somatic variant data for MAVE annotation. Formatting to .vcf"
  ) 

  parser.add_argument(
    "--somatic_variants", 
    type=Path, 
    default="output/variants_with_germline_proximity.tsv",
    help="Path to the variant input file (e.g. output/variants_with_germline_proximity.tsv)"
  )

  parser.add_argument(
     "--neutral_clinvar",
     type=Path, 
     default= "output/neutral_clinvar_filtered.tsv",
     help="Path to the neutral variant file from ClinVar"
  )

  parser.add_argument(
    "--somatic_output",
    type=Path,
    default="output/somatic_variants_GRCh37.vcf",
    help="Path to the annotated output variant file (e.g. output/somatic_variants_GRCh37.vcf)"
  )

  parser.add_argument(
     "--clinvar_output",
     type=Path, 
     default= "output/neutral_clinvar_GRCh37.vcf",
     help="Path to the Clinvar vcf output file."
  )

  return parser.parse_args() 


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

  # Finn ut hva vi skal bruke som ID
  if 'ID' in mapping and mapping['ID'] in df.columns:
      ids = df[mapping['ID']].astype(str)
  else:
      ids = "."

  # create DF 
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
  
  # add EXTRA_INFO if accessible (ClinicalSignificance, ReviewStatus etc.)
  if 'EXTRA_INFO' in mapping:
      for key, col in mapping['EXTRA_INFO'].items():
          vcf_df["INFO"] += f";{key}=" + df[col].astype(str).str.replace(" ", "_")

  # wash data 
  vcf_df = vcf_df.dropna(subset=["REF", "ALT"])
  vcf_df = vcf_df[(vcf_df["REF"] != "-") & (vcf_df["ALT"] != "-")]
  vcf_df = vcf_df[vcf_df["REF"] != vcf_df["ALT"]]
  vcf_df = vcf_df.sort_values(by=["#CHROM", "POS"], ascending=True)
  vcf_df = vcf_df[(vcf_df["REF"].str.len() == 1) & (vcf_df["ALT"].str.len() == 1)]    # keep only SNVs 
  
  final_count = len(vcf_df)
  variants_removed = initial_count - final_count

  print(f"Total variants in: {initial_count}")
  print(f"Variants removed: {variants_removed}")
  print(f"Final variant count: {final_count}\n")

  # write file with header
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

def main():
  
  args = getargs() 

  # SOMATIC VARIANTS 
  somatic_df = pd.read_csv(args.somatic_variants, sep="\t", low_memory=False)
  somatic_mapping = {
    "CHROM": "Chromosome",
    "POS": "Start_Position",
    "ID": ".",
    "REF": "Reference_Allele",
    "ALT": "Tumor_Seq_Allele2",
    "GENE": "Hugo_Symbol"}
  
  tsv_to_vcf(somatic_df, somatic_mapping, args.somatic_output)

  # NEUTRAL CLINVAR VARIANTS 
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