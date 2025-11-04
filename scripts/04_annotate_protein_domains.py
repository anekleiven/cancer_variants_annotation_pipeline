"""
Script: 04_annotate_protein_domains.py
Author: Ane Kleiven

Purpose:
    Annotate cancer variants with protein domain and clan information, 
    and determine whether each variant lies inside or outside a protein domain 

Input:
    - Variants annotated with UniProt accession numbers 
      (e.g., annotation_pipeline/data/variants_with_uniprot.tsv)
    - Pfam domain locations file (Pfam-A.regions.tsv.gz)
    - Pfam clan metadata file (Pfam-A.clans.tsv)

Output:
    - TSV file of cancer variants annotated with Pfam domains and clans
      (e.g., annotation_pipeline/output/variants_with_pfam_domains.tsv)
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
        required=True,
        help="Path to the input file (e.g. annotation_pipeline/data/variants_with_uniprot.tsv)"
    )
    
    parser.add_argument(
        "--pfam",
        type=Path,
        required=True,
        help="Path to the Pfam domain locations file (e.g. annotation_pipeline/data/Pfam-A.regions.tsv.gz)"
    )
    
    parser.add_argument(
        "--pclan",
        type=Path,
        required=True,
        help="Path to the Pfam clan file (e.g. annotation_pipeline/data/Pfam-A.clans.tsv)"
    )
    
    parser.add_argument(
        "--output",
        type=Path, 
        required=False,
        default="annotation_pipeline/output/variants_with_pfam_domains.tsv",
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

def load_pfam_regions(pfam_regions_file: Path, chunksize: int = 500000) -> pd.DataFrame:
    """Load Pfam-A regions file efficiently and clean duplicates."""
    print(f"Loading Pfam-A domain regions: {pfam_regions_file}")

    use_cols = ["pfamseq_acc", "pfamA_acc", "seq_start", "seq_end"]
    chunks = []

    for chunk in pd.read_csv(pfam_regions_file, sep="\t", 
                             usecols=use_cols, 
                             chunksize=chunksize):
        chunk = chunk.drop_duplicates(subset=use_cols)
        chunks.append(chunk)

    df = pd.concat(chunks).drop_duplicates(subset=use_cols)

    df = df.rename(columns={
        "pfamseq_acc": "UNIPROT_ACCESSION",
        "pfamA_acc": "PFAM_ACCESSION",
        "seq_start": "DOMAIN_START",
        "seq_end": "DOMAIN_END"
        }
        )

    df["DOMAIN_START"] = pd.to_numeric(df["DOMAIN_START"], errors="coerce").astype("Int64")
    df["DOMAIN_END"] = pd.to_numeric(df["DOMAIN_END"], errors="coerce").astype("Int64")

    print(f"Loaded {len(df):,} unique Pfam domain mappings.\n")
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
    """
    Merge cancer variants with Pfam domain and clan metadata,
    and determine whether each variant lies within a protein domain.

    Output: one row per variant, with overlapping domain(s) or 'NA'.
    """
    
    print(f"Merging variants with Pfam domain regions...\n")
    
    print("Step 1: Merge variants with Pfam regions")
    merged = variants.merge(pfam_regions, on="UNIPROT_ACCESSION", how="left")
    print(f"After merge: {len(merged):,} rows.\n")

    print("Step 2: Add Pfam domain names and clans")
    merged = merged.merge(pfam_clans, on="PFAM_ACCESSION", how="left")
    print(f"After adding pfam clan information: {len(merged):,} rows.\n")

    print("Step 3: Determine whether variant lie within domains")

    # Check data types (numeric columns) 
    for col in ["Protein_position", "DOMAIN_START", "DOMAIN_END"]:
        merged[col] = pd.to_numeric(merged[col], errors='coerce')
    
    # Check if variant position overlaps with protein domain 
    merged["IN_DOMAIN"] = (
        (merged["Protein_position"] >= merged['DOMAIN_START']) & 
        (merged["Protein_position"] <= merged['DOMAIN_END'])
    )

    print("Example overlap check:\n",
          merged[["Hugo_Symbol", "HGVSp_Short", "Protein_position",
                  "DOMAIN_NAME", "DOMAIN_START", "DOMAIN_END", 
                  "IN_DOMAIN"]]
          .head(5), "\n")
    
    print("Step 4: Summarize per variant")
    # Only keep overlapping domain rows 
    overlapping = merged[merged["IN_DOMAIN"]].copy() 

    # Collapse multiple overlapping domains if they exist 
    # Groups all rows that belong to the same variant (gene + protein change) 
    # Applies an aggregation function that removes duplicated domain names, sorts them alphabetical and joins them to one string
    summarized = (
        overlapping.groupby(["Hugo_Symbol", "HGVSp_Short"], as_index=False)
        .agg({
            "DOMAIN_NAME": lambda x: ";".join(sorted(set(x.dropna()))),
            "DESCRIPTION": lambda x: ";".join(sorted(set(x.dropna())))
            })
            )
    
    # Keep original variant order, fill missing domains with NA 
    out = (
        variants.merge(
            summarized, 
            on=["Hugo_Symbol", "HGVSp_Short"], 
            how="left"))
    
    out["DOMAIN_NAME"] = out["DOMAIN_NAME"].fillna("NA")


    print("Step 5: Final annotated dataset")
    print(out.head(), "\n")
    print(f"Output: {len(out):,} variants with protein domain annotations.\n")

    return out

# ============================================================
# 5. Main function 
# ============================================================

def main():
    args = get_args() 

    variants = load_variants_with_uniprot(args.input)
    pfam_domains = load_pfam_regions(args.pfam)
    pfam_clans = load_pfam_clans(args.pclan)

    annotated = merge_variants_with_pfam(variants, pfam_domains, pfam_clans)

    annotated.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved annotated variant file to: {args.output}\n")


# ============================================================
# 6. Entry point
# ============================================================

if __name__ == "__main__":
    main()
