""" 
====================================================================
Somatic Variant file (.tsv) to VCF format 
Preparation for MAVE annotation 
====================================================================

Script: 06_mave_tsv_to_vcf.py 
Author: Ane Kleiven 

Major outputs: 
  Prepare somatic variant data for MAVE annotation using VEP
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
    "--input", 
    type=Path, 
    default="output/variants_with_germline_proximity.tsv",
    help="Path to the variant input file (e.g. output/variants_with_germline_proximity.tsv)"
  )

  parser.add_argument(
    "--output",
    type=Path,
    default="output/somatic_variants_GRCh37.vcf",
    help="Path to the annotated output variant file (e.g. output/somatic_variants_GRCh37.vcf)"
  )

  return parser.parse_args() 


def main():

  args = getargs() 

  # load somatic variant file
  print("\nLoading somatic variant data...\n")
  somatic_variants = pd.read_csv(args.input, sep="\t", low_memory=False)
  print(f"Loaded {len(somatic_variants)} variants from {args.input.name}\n")

  # convert file to VCF format for MAVE VEP annotation 
  print("Converting file to VCF format for MAVE annotation...\n")
  vcf_df = pd.DataFrame({
    "#CHROM": "chr" + somatic_variants["Chromosome"].astype(str),
    "POS": somatic_variants["Start_Position"].astype(int),
    "ID": ".",
    "REF": somatic_variants["Reference_Allele"],
    "ALT": somatic_variants["Tumor_Seq_Allele2"],
    "QUAL": ".",
    "FILTER": somatic_variants["FILTER"].fillna("PASS"),
    "INFO": (
        "GENE=" + somatic_variants["Hugo_Symbol"].astype(str) +
        ";ENTREZ=" + somatic_variants["Entrez_Gene_Id"].fillna(0).astype(int).astype(str)
    )
})

  # Remove rows where reference allele or alternate allele is missing 
  print("Removing rows where reference allele or alternate allele is missing...\n")
  vcf_df = vcf_df.dropna(subset=["REF", "ALT"])
  print("Removing rows where reference allele or alternate allelse is a dash ('-')")
  vcf_df = vcf_df[(vcf_df["REF"] != "-") & (vcf_df["ALT"] != "-")]
  
  # Remove rows where reference allele is equal to alternate allele 
  print("Removing rows where reference allele is equal to alternate allele...\n")
  vcf_df = vcf_df[vcf_df["REF"] != vcf_df["ALT"]]

  # Add VCF header 
  print("Adding VCF header...\n")
  vcf_header = [
      "##fileformat=VCFv4.2",
      "##source=Master_thesis_anekleiven",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  ]

  # sort by Chromosome and Location (ascending)
  print("Sorting variants by Chromosome and Position...\n")
  vcf_df = vcf_df.sort_values(by=["#CHROM", "POS"], ascending=True)

  print("Preview of the VCF variant file:")
  print(vcf_df.head(),"\n")

  # Save vcf file to output folder 
  output_vcf = args.output

  print("Writing file and saving to output folder...\n")
  with open(output_vcf, "w") as f:
      for line in vcf_header:
          f.write(line + "\n")
      vcf_df.to_csv(f, sep="\t", index=False, header=False)
  print("VCF saved as output/somatic_variants_GRCh37.vcf")

  print("Script ended!🎉\n")

if __name__ == "__main__":
   main() 