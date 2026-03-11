"""
====================================================================
Protein Domain Annotation Script
====================================================================

Script: 03_annotate_protein_domains.py
Author: Ane Kleiven

Purpose:
    Annotate cancer variants with protein domain and clan information, 
    and determine whether each variant lies inside or outside a protein domain 

Input:
    - Variants annotated with UniProt accession numbers 
      (e.g., data/variants_with_uniprot.tsv)
    - Pfam domain locations file (Pfam-A.regions.tsv.gz)
    - Pfam clan metadata file (Pfam-A.clans.tsv)

Output:
    - TSV file of cancer variants annotated with Pfam domains and clans
      (e.g., output/variants_with_pfam_domains.tsv)
"""

import pandas as pd
import argparse
from pathlib import Path


# ============================================================
# Parse command line arguments 
# ============================================================

def get_args(): 
    parser = argparse.ArgumentParser(
        description="Annotate cancer variants with protein domains and clan information"
    )

    parser.add_argument(
        "--input",
        type=Path,
        required=False,
        default="data/variants_with_uniprot.tsv",
        help="Path to the input file (e.g. data/variants_with_uniprot.tsv)"
    )
    
    parser.add_argument(
        "--pfam",
        type=Path,
        required=False,
        default="data/Pfam-A.regions.tsv.gz",
        help="Path to the Pfam domain locations file (e.g. data/Pfam-A.regions.tsv.gz)"
    )
    
    parser.add_argument(
        "--pclan",
        type=Path,
        required=False,
        default="data/Pfam-A.clans.tsv",
        help="Path to the Pfam clan file (e.g. data/Pfam-A.clans.tsv)"
    )
    
    parser.add_argument(
        "--output",
        type=Path, 
        required=False,
        default="output/variants_with_pfam_domains.tsv",
        help="Path to save the annotated output TSV file"
    )
    
    return parser.parse_args()


# ============================================================
# 1. Load variants with UniProt accessions
# ============================================================

def load_variants_with_uniprot(input_file: Path) -> pd.DataFrame:
    """Load variant table and standardize column names and data types."""
    print(f"\nLoading variant file: {input_file}")

    df = pd.read_csv(input_file, sep="\t", low_memory=False)
    df = df.rename(columns={"uniprot_acc": "UNIPROT_ACCESSION"})

    int_cols = [
        "Entrez_Gene_Id", "Start_Position", "End_Position",
        "t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count",
        "n_depth", "t_depth", "Protein_position", "Exon_Number",
        "Score", "Samples"
    ]

    for col in int_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    print(f"Loaded {len(df):,} variant rows.\n")
    return df


# ============================================================
# 2. Load Pfam domain locations (Pfam-A.regions.tsv.gz)
# ============================================================

def load_pfam_regions(pfam_regions_file: Path, relevant_uniprot: list, chunksize: int = 500000) -> pd.DataFrame:
    """Load Pfam-A regions file and filter for relevant proteins during chunking."""
    print(f"Loading and filtering Pfam-A domain regions: {pfam_regions_file}")

    use_cols = ["pfamseq_acc", "pfamA_acc", "seq_start", "seq_end"]
    chunks = []

    # Reading in chunks to prevent memory to explode 
    for chunk in pd.read_csv(pfam_regions_file, sep="\t", 
                             usecols=use_cols, 
                             chunksize=chunksize):
        
        # Drop rows with other pfam accession numbers (not in our data)
        chunk = chunk[chunk["pfamseq_acc"].isin(relevant_uniprot)]
        
        if not chunk.empty:
            chunk = chunk.drop_duplicates(subset=use_cols)
            chunks.append(chunk)

    if not chunks:
        print("No matching UniProt accessions found in Pfam file.")
        return pd.DataFrame(columns=["UNIPROT_ACCESSION", "PFAM_ACCESSION", "DOMAIN_START", "DOMAIN_END"])

    df = pd.concat(chunks).drop_duplicates(subset=use_cols)

    df = df.rename(columns={
        "pfamseq_acc": "UNIPROT_ACCESSION",
        "pfamA_acc": "PFAM_ACCESSION",
        "seq_start": "DOMAIN_START",
        "seq_end": "DOMAIN_END"
        })

    # Convert to numeric data type 
    df["DOMAIN_START"] = pd.to_numeric(df["DOMAIN_START"], errors="coerce").astype("Int64")
    df["DOMAIN_END"] = pd.to_numeric(df["DOMAIN_END"], errors="coerce").astype("Int64")

    print(f"Filtered Pfam down to {len(df):,} relevant domain mappings.\n")
    return df

# ============================================================
# 3. Load Pfam clans 
# ============================================================

def load_pfam_clans(pfam_clans_file: Path) -> pd.DataFrame:
    """Load Pfam clan metadata."""
    print(f"Loading Pfam clans: {pfam_clans_file}")
    df = pd.read_csv(
        pfam_clans_file, 
        sep="\t", 
        header=None,
        names=["PFAM_ACCESSION", "CLAN_ACCESSION", "CLAN_ID", "DOMAIN_NAME", "DESCRIPTION"]
    )

    print(f"Loaded {len(df):,} Pfam clan entries.\n")
    return df


# ============================================================
# 4. Merge everything
# ============================================================

def merge_variants_with_pfam(
        variants: pd.DataFrame, 
        pfam_regions: pd.DataFrame, 
        pfam_clans: pd.DataFrame
) -> pd.DataFrame:
    
    print(f"Merging variants with Pfam domain regions...\n")
    
    # 1. Get unique accessions present in your variants to shrink the Pfam table
    # This is the most important step for saving RAM
    relevant_uniprot = variants["UNIPROT_ACCESSION"].unique()
    pfam_regions = pfam_regions[pfam_regions["UNIPROT_ACCESSION"].isin(relevant_uniprot)].copy()
    
    print(f"Filtered Pfam regions down to {len(pfam_regions):,} relevant records.")

    # 2. Perform the merge
    merged = variants.merge(pfam_regions, on="UNIPROT_ACCESSION", how="inner")
    
    # 3. Filter for overlap 
    # Using .query() for memory efficiency 
    overlapping = merged[
        (merged["Protein_position"] >= merged['DOMAIN_START']) & 
        (merged["Protein_position"] <= merged['DOMAIN_END'])
    ].copy()
    
    # Clean up the massive 'merged' object from RAM 
    del merged 

    print(f"Found {len(overlapping):,} overlapping domain instances.")

    # 4. Add Clan info to the overlapping rows
    overlapping = overlapping.merge(pfam_clans, on="PFAM_ACCESSION", how="left")

    # 5. Summarize and collapse multiple domains per variant 
    summarized = (
        overlapping.groupby(["Hugo_Symbol", "HGVSp_Short"], as_index=False)
        .agg({
            "DOMAIN_NAME": lambda x: ";".join(sorted(set(x.dropna()))),
            "DESCRIPTION": lambda x: ";".join(sorted(set(x.dropna())))
        })
    )
    
    # 6. Join back to the original variant list to restore 'NA' rows
    out = variants.merge(summarized, on=["Hugo_Symbol", "HGVSp_Short"], how="left")
    out["DOMAIN_NAME"] = out["DOMAIN_NAME"].fillna("NA")

    # 7. Add IN_DOMAIN binary column 
    out["IN_DOMAIN"] = (out["DOMAIN_NAME"] != "NA") & (out["DOMAIN_NAME"] != "")

    return out

# ============================================================
# 5. Main function 
# ============================================================

def main():
    args = get_args() 

    # 1. Load variants first
    variants = load_variants_with_uniprot(args.input)
    
    # 2. Get list of UniProt IDs to use as a filter
    relevant_uniprot = variants["UNIPROT_ACCESSION"].unique().tolist()

    # 3. Load Pfam using filter
    pfam_domains = load_pfam_regions(args.pfam, relevant_uniprot)
    
    pfam_clans = load_pfam_clans(args.pclan)

    annotated = merge_variants_with_pfam(variants, pfam_domains, pfam_clans)

    annotated.to_csv(args.output, sep="\t", index=False)

    print(f"\nSaved annotated variant file to: {args.output}\n")

# ============================================================
# 6. Entry point
# ============================================================

if __name__ == "__main__":
    main()
