"""
Script: 02_annotate_hotspots.py
Author: Ane Kleiven

Script to annotate cancer hotspots to 
cancer variants

Input: 
	OncoKB annotated cancer variants (.MAF file) 
	Hotspots file, long format (e.g. hotspots_long.tsv)

Join tables by: 
'Hugo_Symbol' and 'HGVSp_Short'
"""

import pandas as pd
import argparse 
from pathlib import Path

# define function for user input file paths 
def get_args():
	parser = argparse.ArgumentParser(description="Annotate variants with cancer hotspots")
	parser.add_argument(
		"--input",
		type=Path,
		required=True, 
		help="Path to the input file (e.g.annotated_output.maf)"
  )
	
	parser.add_argument(
    "--hotspots",
      type=Path,
      required=True,
      help="Path to the expanded cancer hotspots file"
	)

	parser.add_argument(
		"--output",
		type=Path,
		required=False,
		default=Path("annotation_pipeline/output/annotated_with_hotspots.maf"),
		help="Path to save the annotated output file (default: annotation_pipeline/output/annotated_with_hotspots.maf)"
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
	
	# merge by gene and protein change: 
	merged = pd.merge(
		variants, 
		hotspots,
		how="left",
		left_on=['Hugo_Symbol', 'HGVSp_Short'],
		right_on=['Hugo_Symbol', 'HGVSp_Short'],
	)

	print(f"Merged dataset has {len(merged):,} rows\n") 
	
	# save to output 
	merged.to_csv(args.output, sep="\t", index=False)
	print(f"Annotated file saved to:\n {args.output.resolve()}🥳\n")


if __name__ == "__main__":
	main() 

