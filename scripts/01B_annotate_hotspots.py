
# ====================================================================
# Cancer Hotspots Annotation 
# ====================================================================

"""
Script: 01B_annotate_hotspots.py
Author: Ane Kleiven

Script to annotate cancer hotspots to cancer variants

Input: 
	OncoKB annotated cancer variants (.MAF file) 
	Hotspots variant file (long format: one row per variant), see '01A_hotspots_long.py' 

Join tables by: 
'Hugo_Symbol' and 'HGVSp_Short'

"""

import pandas as pd
import argparse 
from pathlib import Path

def get_args():
	parser = argparse.ArgumentParser(
		description="Annotate variants with cancer hotspots"
		)

	parser.add_argument(
		"--input",
		type=Path,
		required=True, 
		help="Path to the input file (e.g. annotated_output.maf)"
	)
	
	parser.add_argument(
		"--hotspots",
		type=Path,
		required=True, 
		help="Path to the cancer hotspots file (e.g. hotspots_long.tsv)"
	)

	parser.add_argument(
		"--output",
		type=Path,
		required=True,
		help="Path to save the annotated output file (e.g. annotated_with_hotspots.tsv)."
	)
	
	return parser.parse_args() 


def main():
	# Call get_args() to get the users input file path: 
	args=get_args()

	# Load files 
	print(f"\nLoading files..\n")
	variants = pd.read_csv(args.input, sep="\t", low_memory=False) 
	hotspots = pd.read_csv(args.hotspots, sep="\t", low_memory=False)

	print(f"Loaded {len(variants):,} variants from {args.input.name}") 
	print(f"Loaded {len(hotspots):,} hotspots from {args.hotspots.name}\n")

	# Add In_Hotspot boolean column
	hotspots["In_Hotspot"] = True

	# Merge hotspots file with variants file
	merged = pd.merge(
		variants, 
		hotspots,
		how="left",
		on=["Hugo_Symbol", "HGVSp_Short"]
	)

	merged["In_Hotspot"] = merged["In_Hotspot"].fillna(False) 

	print(f"Variants in hotspot: {merged['In_Hotspot'].sum():,}")

	# Save output
	merged.to_csv(args.output, sep="\t", index=False)
	print(f"\nAnnotated file saved to:\n {args.output.resolve()} 🥳\n")


	print("Cancer hotspots annotation complete!🥳")


if __name__ == "__main__":
	main()