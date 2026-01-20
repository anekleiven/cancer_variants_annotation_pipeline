"""
====================================================================
Hotspot file formating script
====================================================================

Script: 01_hotspots_long.py
Author: Ane Kleiven

Script to transform the hotspots.txt file into a long format, 
separating the different protein changes.

hotspots.txt file can be downloaded here: 
https://www.cancerhotspots.org/#/home

Example usage: 
    python hotspots_long.py \
        --input data/hotspots.txt \
        --output output/hotspots_long.tsv
"""

import pandas as pd 
import argparse
from pathlib import Path

def get_args(): 
    parser= argparse.ArgumentParser(
        description="Transform cancer hotspots file into long format"
    ) 

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to the input hotspot file (e.g. data/hotspots.txt)"
    )

    parser.add_argument(
        "--output",
        type=Path, 
        required=False,
        default="data/hotspots_long.tsv",
        help="Path to the output hotspot file (defalt: data/hotspots_long.tsv)"
    )

    return parser.parse_args()


def main():
    args=get_args() 

    print(f"\nTransforming hotspots file into long format:\n")
    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}\n")


    # Load the input file
    print(f"Reading hotspots file...\n")
    hotspots = pd.read_csv(args.input, sep="\t", header=0)

    expanded_rows = [] 

    # Iterate over each row in the original file
    print(f"Splitting variants into separate rows.\n")
    for _, row in hotspots.iterrows():
        gene = row['Gene'] 
        residue = row['Residue']
        hotspot_type = row['Type']          # renamed 'type' (since type() is a Python keyword)
        qvalue = row['Q-value']
        samples = row['Samples'] 
        tumor_comp = row['Tumor Type Composition']   # FIXED: previously you had [‘Tumor Type Composition’] with brackets!

        # Split the variants like "R:204|K:142|L:46"
        for variant in str(row['Variants']).split('|'):
            alt_aa = variant.split(':')[0]  # get the amino acid letter before the colon
            protein_change = f"p.{residue}{alt_aa}"

            expanded_rows.append({
                'Hugo_Symbol': gene,
                'HGVSp_Short': protein_change,
                'Hotspot_Type': hotspot_type,
                'Q_Value': qvalue,
                'Samples': samples, 
                'Tumor_Types': tumor_comp
            })

    # Convert to DataFrame
    hotspots_long = pd.DataFrame(expanded_rows)
    print(f"File expanded to {len(hotspots_long):,} rows\n")

    # Save as a .txt file 
    hotspots_long.to_csv(args.output, sep="\t",index=False)
    print(f"File saved to: {args.output.resolve()}")

    # Preview output file 
    print(f"\nPreview of the first ten rows:\n")
    print(f"{hotspots_long.head(10)}\n")


if __name__ == "__main__":
    main() 