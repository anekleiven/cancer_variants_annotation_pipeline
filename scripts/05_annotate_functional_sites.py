""" 

Script: 05_annotate_functional_sites.py 
Author: Ane Kleiven 

Uniprot API: 
    For each protein in the variant data set:
    Ask Uniprot for the full feature annotation (JSON) 
    Pick out relevant functional sites (in ftype) 
    Extract residue positions and descriptions 

Merge cancer variants with functional sites df 
Check if variant is inside or outside a functional site 

"""
import argparse
from pathlib import Path
import requests
import pandas as pd
import time
from tqdm import tqdm

def get_args(): 
    parser = argparse.ArgumentParser(
        description="Annotate cancer variants with protein functional sites"
    )

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to the input variant file (e.g. annotation_pipeline/output/variants_with_pfam_domains.tsv)"
    )

    parser.add_argument(
        "--funcsites",
        type=Path,
        required=False,
        default="annotation_pipeline/data/human_functional_sites.tsv",
        help="Path to the functional sites file (e.g. annotation_pipeline/data/human_functional_sites.tsv)"
    )

    parser.add_argument(
        "--output",
        type=Path, 
        required=False,
        default="annotation_pipeline/output/variants_with_func_sites.tsv",
        help="Path to the output variant file (e.g. annotation_pipeline/output/variants_with_func_sites.tsv)"
    )

    return parser.parse_args()


def get_functional_sites(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    r = requests.get(url)
    if r.status_code != 200:
        return []
    data = r.json()
    features = []
    for f in data.get("features", []):
        ftype = f.get("type")
        if ftype in [
            "Active site",
            "Binding site",
            "Metal binding", 
            "Modified residue", 
            "Region",
            "Motif",
            "Site",
            "Signal peptide", 
            "Transmembrane",
            "Topological domain",
            "Glycosylation site", 
        ]:
            loc = f.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value", start)    # sets end to start if end is missing 
            desc = f.get("description", "")
            features.append({                               # stores each extracted feature in a dictionary 
                "UNIPROT_ACCESSION": uniprot_id,
                "FEATURE_TYPE": ftype,
                "FUNC_SITE_DESC": desc,
                "FUNC_SITE_START": start,
                "FUNC_SITE_END": end
            })
    return features


def merge_functional_sites(variants, functional_sites_df) -> pd.DataFrame:
    """Merge variants with functional sites downloaded from Uniprot"""
    print(f"Merging variants with Uniprot functional sites...\n ")
    merged = variants.merge(functional_sites_df,on="UNIPROT_ACCESSION", how="left")
    
    print(f"Final annotated dataset: {len(merged):,} rows.\n")
    return merged


def in_functional_site(df: pd.DataFrame) -> pd.DataFrame:
    """
    Determine if each variant lies inside or outside a functional site.
    Adds:
        - IN_FUNC_SITE (True/False)
        - FEATURE_TYPE (semicolon-joined site types)
        - FUNC_SITE_DESC (semicolon-joined site descriptions)
    """
    print("Determining if variants are inside or outside functional site...\n")

    # Ensure numeric positions
    df["Protein_position"] = pd.to_numeric(df["Protein_position"], errors="coerce")
    df["FUNC_SITE_START"] = pd.to_numeric(df["FUNC_SITE_START"], errors="coerce")
    df["FUNC_SITE_END"] = pd.to_numeric(df["FUNC_SITE_END"], errors="coerce")

    # Boolean overlap
    df["IN_FUNC_SITE_tmp"] = (
        (df["Protein_position"] >= df["FUNC_SITE_START"]) &
        (df["Protein_position"] <= df["FUNC_SITE_END"])
    )

    # Summarize only overlapping variants (far fewer rows)
    overlapping = df[df["IN_FUNC_SITE_tmp"]].copy()
    summarized = (
        overlapping.groupby(["Hugo_Symbol", "HGVSp_Short"], as_index=False)
        .agg({
            "FEATURE_TYPE": lambda x: ";".join(sorted(set(x.dropna()))),
            "FUNC_SITE_DESC": lambda x: ";".join(sorted(set(x.dropna())))
        })
    )

    # Merge with the unique variant-level table (one row per variant)
    base_variants = (
        df[["Hugo_Symbol", "HGVSp_Short"] + [c for c in df.columns if c not in [
            "Hugo_Symbol", "HGVSp_Short", "FEATURE_TYPE", "FUNC_SITE_DESC",
            "FUNC_SITE_START", "FUNC_SITE_END", "IN_FUNC_SITE_tmp"
        ]]]
        .drop_duplicates(subset=["Hugo_Symbol", "HGVSp_Short"])
    )

    annotated = base_variants.merge(
        summarized,
        on=["Hugo_Symbol", "HGVSp_Short"],
        how="left"
    )

    annotated["FEATURE_TYPE"] = annotated["FEATURE_TYPE"].fillna("NA")
    annotated["FUNC_SITE_DESC"] = annotated["FUNC_SITE_DESC"].fillna("NA")
    annotated["IN_FUNC_SITE"] = annotated["FEATURE_TYPE"].ne("NA")

    print("Example summarized annotations:\n")
    print(annotated[["Hugo_Symbol", "HGVSp_Short", "IN_FUNC_SITE",
                     "FEATURE_TYPE", "FUNC_SITE_DESC"]].head(), "\n")

    return annotated



def main(): 
    args = get_args() 

    # Load variants
    print("\nLoading variants file...\n")
    variants = pd.read_csv(args.input, sep="\t", low_memory=False)

    # Fetch or load functional sites
    accessions = variants["UNIPROT_ACCESSION"].dropna().unique()

    print("Searching for functional sites file...\n")
    if Path(args.funcsites).exists(): 
        print(f"Functional sites file found: {args.funcsites}, loading instead of downloading.\n")
        functional_sites_df = pd.read_csv(args.funcsites, sep="\t", low_memory=False)
    else: 
        print("Functional sites file not found.")
        print("Fetching functional sites from UniProt API...\n")
        functional_sites = []
        for acc in tqdm(accessions, desc="Fetching functional sites"):
            functional_sites.extend(get_functional_sites(acc))     # extend() appends each proteins features to one big list 
            time.sleep(0.25)                                       # avoids overloading Uniprots server      
        functional_sites_df = pd.DataFrame(functional_sites)
        functional_sites_df.to_csv(args.funcsites, sep="\t", index=False)

    for col in ["FUNC_SITE_START", "FUNC_SITE_END"]:
        functional_sites_df[col] = pd.to_numeric(functional_sites_df[col], errors="coerce").astype("Int64")

    print(f"Example output from 'functional_sites_df':")
    print(f"\n{functional_sites_df.head()}")

    # Merge + annotate 
    merged_variants_df = merge_functional_sites(variants, functional_sites_df) 
    annotated = in_functional_site(merged_variants_df)

    # Summary 
    n_total = len(annotated)
    n_true = annotated["IN_FUNC_SITE"].sum()
    print(f"{n_true:,} of {n_total:,} variants overlap functional sites ({n_true/n_total:.1%}).\n")

    # Save final annotated file 
    annotated.to_csv(args.output, sep="\t", index=False) 
    print(f"\nAnnotated variants file saved to: {args.output}\n")


if __name__ == "__main__":
    main()

