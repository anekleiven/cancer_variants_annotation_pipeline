"""
====================================================================
Cancer Hotspots Annotation 
====================================================================

Script: 01_annotate_hotspots.py
Author: Ane Kleiven

Script to annotate cancer hotspots to cancer variants

Input: 
	OncoKB annotated cancer variants (.MAF file) 
	Hotspots variant file (V2 from 2017) 

Join tables by: 
'Hugo_Symbol' and 'HGVSp_Short'
'Hugo_Symbol' and 'Codon' 
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
		"--hotspots_v2",
		type=Path,
		default="/home/anekl/git/master/cancer_variants_annotation_pipeline/data/cancerhotspots.v2.maf",
		required=False, 
		help="Path to the v2 cancer hotspots file."
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
	args = get_args()

	# load variant file  
	print(f"\nLoading variant file..🤓")
	variants = pd.read_csv(args.input, sep="\t", low_memory=False) 
	print(f"Loaded {len(variants):,} variants from {args.input.name}") 

	# extract codon from HGVSp_Short before any merging
	variants['Codon'] = variants["HGVSp_Short"].str.extract(r'p\.([A-Z]\d+)')

	# get unique gene + protein change combinations for filtering V2
	variant_keys = set(zip(variants['Hugo_Symbol'], variants['HGVSp_Short']))

	# --------------------------------------------------------
	# Load V2 hotspots
	# --------------------------------------------------------

	print(f"\nLoading V2 hotspots file..🤓")
	chunks = [] 

	for chunk in pd.read_csv(args.hotspots_v2, sep="\t", comment="#", low_memory=False, chunksize=50000):
		keep = [
			(gene, hgvsp) in variant_keys
			for gene, hgvsp in zip(chunk['Hugo_Symbol'], chunk['HGVSp_Short'])
		]
		chunks.append(chunk[keep])

	hotspots_v2_filtered = pd.concat(chunks).drop_duplicates(subset=['Hugo_Symbol', 'HGVSp_Short'])
	hotspots_v2_filtered.to_csv("data/hotspots_v2_filtered.tsv", sep="\t", index=False)
	print(f"Filtered to {len(hotspots_v2_filtered):,} unique V2 hotspot variants")

	hotspots_v2_filtered["is_hotspot"] = True

	# merge V2 hotspots with variants file
	merged = pd.merge(
		variants, 
		hotspots_v2_filtered[["Hugo_Symbol", "HGVSp_Short", "is_hotspot"]],
		how="left",
		on=["Hugo_Symbol", "HGVSp_Short"]
	)
	merged["is_hotspot"] = merged["is_hotspot"].fillna(False)


	# verify no new rows were added
	assert len(merged) == len(variants), f"Row count changed! {len(variants):,} -> {len(merged):,}"
	print(f"\nRow count verified: {len(merged):,} variants")

	# summary
	print(f"Variants in hotspot: {merged['is_hotspot'].sum():,}")

	# save output
	merged.to_csv(args.output, sep="\t", index=False)
	print(f"\nAnnotated file saved to:\n {args.output.resolve()} 🥳\n")


if __name__ == "__main__":
	main()