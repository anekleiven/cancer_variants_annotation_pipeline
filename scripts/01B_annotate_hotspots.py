"""
====================================================================
Cancer Hotspots Annotation 
====================================================================

Script: 01B_annotate_hotspots.py
Author: Ane Kleiven

Script to annotate cancer hotspots to cancer variants

Input: 
	OncoKB annotated cancer variants (.MAF file) 
	Hotspots variant file (long format: one row per variant)

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
		default="/home/anekl/git/master/genie_oncokb_processing_scripts/data/annotated_output.maf",
		required=False, 
		help="Path to the input file (e.g.annotated_output.maf)"
	)
	
	parser.add_argument(
		"--hotspots",
		type=Path,
		default="/home/anekl/git/master/cancer_variants_annotation_pipeline/data/hotspots_long.tsv",
		required=False, 
		help="Path to the cancer hotspots file (long format)"
	)

	parser.add_argument(
		"--output",
		type=Path,
		required=False,
		default=Path("/home/anekl/git/master/cancer_variants_annotation_pipeline/output/annotated_with_hotspots.maf"),
		help="Path to save the annotated output file."
	)
	
	return parser.parse_args() 


def main():
	# call get_args() to get the users input file path: 
	args=get_args()

	# load files 
	print(f"\nLoading files...🤓\n")
	variants = pd.read_csv(args.input, sep="\t", low_memory=False) 
	hotspots = pd.read_csv(args.hotspots, sep="\t", low_memory=False)

	print(f"Loaded {len(variants):,} variants from {args.input.name}") 
	print(f"Loaded {len(hotspots):,} hotspots from {args.hotspots.name}\n")

	# add in_hotspot boolean 
	hotspots["in_hotspot"] = True

	# merge hotspots file with variants file
	merged = pd.merge(
		variants, 
		hotspots,
		how="left",
		on=["Hugo_Symbol", "HGVSp_Short"]
	)

	merged["in_hotspot"] = merged["in_hotspot"].fillna(False) 


	# summary
	print(f"Variants in hotspot: {merged['in_hotspot'].sum():,}")

	# save output
	merged.to_csv(args.output, sep="\t", index=False)
	print(f"\nAnnotated file saved to:\n {args.output.resolve()} 🥳\n")


	print("Cancer hotspots annotation complete!🥳")


if __name__ == "__main__":
	main()