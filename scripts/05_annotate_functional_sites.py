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
            end = loc.get("end", {}).get("value", start)
            desc = f.get("description", "")
            features.append({
                "UNIPROT_ACCESSION": uniprot_id,
                "FEATURE_TYPE": ftype,
                "FUNC_SITE_DESC": desc,
                "FUNC_SITE_START": start,
                "FUNC_SITE_END": end
            })
    return features


def find_overlapping_sites(row, sites_dict):
    """
    For a single variant, find all overlapping functional sites.
    Checking each variant against the dictionary of functional sites. 
    Returns semicolon-separated strings or NA.
    """
    # pull the protein accession and amino acid position for a variant 
    uniprot_acc = row["UNIPROT_ACCESSION"]
    position = row["Protein_position"]
    
    # if either uniprot_acc or position is missing, there can't be an overlap 
    if pd.isna(uniprot_acc) or pd.isna(position):
        return "NA", "NA", False
    
    # Get functional sites for this protein
    sites = sites_dict.get(uniprot_acc, [])
    
    # No sites = no overlap 
    if not sites:
        return "NA", "NA", False
    
    # Collect all overlapping types and description for the variant 
    overlapping_types = []
    overlapping_descs = []
    
    # Iterates over each site and extracts its start/end position 
    for site in sites:
        start = site["FUNC_SITE_START"]
        end = site["FUNC_SITE_END"]
        
        # Check if variants position overlaps 
        if pd.notna(start) and pd.notna(end):
            if start <= position <= end:
                # If overlap, append the site 
                if pd.notna(site["FEATURE_TYPE"]):
                    overlapping_types.append(str(site["FEATURE_TYPE"]))
                if pd.notna(site["FUNC_SITE_DESC"]):
                    overlapping_descs.append(str(site["FUNC_SITE_DESC"]))
    
    # If overlaps: remove duplicates (set), sort alphabetical (sorted), collapse to one string (.join) 
    if overlapping_types:
        types_str = ";".join(sorted(set(overlapping_types)))
        descs_str = ";".join(sorted(set(overlapping_descs))) if overlapping_descs else "NA"
        return types_str, descs_str, True
    # No overlaps:
    else:
        return "NA", "NA", False


def annotate_with_functional_sites(variants, functional_sites_df) -> pd.DataFrame:
    """
    Memory-efficient annotation: iterate through variants and check overlaps.
    """
    print(f"Annotating {len(variants):,} variants with functional sites...\n")
    
    # Ensure numeric types
    variants["Protein_position"] = pd.to_numeric(variants["Protein_position"], errors="coerce")
    functional_sites_df["FUNC_SITE_START"] = pd.to_numeric(functional_sites_df["FUNC_SITE_START"], errors="coerce")
    functional_sites_df["FUNC_SITE_END"] = pd.to_numeric(functional_sites_df["FUNC_SITE_END"], errors="coerce")
    
    # Create a dictionary: {UNIPROT_ACCESSION: [list of sites]}
    # For memory efficiency 
    print("Building functional sites dictionary...\n")
    sites_dict = {}
    for _, row in functional_sites_df.iterrows():
        acc = row["UNIPROT_ACCESSION"]
        if acc not in sites_dict:
            sites_dict[acc] = []
        sites_dict[acc].append({
            "FEATURE_TYPE": row["FEATURE_TYPE"],
            "FUNC_SITE_DESC": row["FUNC_SITE_DESC"],
            "FUNC_SITE_START": row["FUNC_SITE_START"],
            "FUNC_SITE_END": row["FUNC_SITE_END"]
        })
    
    # Find the number of accessions with functional sites
    print(f"Built lookup dictionary for {len(sites_dict):,} proteins.\n")
    
    # Apply the 'find_overlapping_sites' to every variant 
    print("Finding overlapping sites for each variant...\n")
    results = variants.apply(
        lambda row: find_overlapping_sites(row, sites_dict),
        axis=1
    )
    
    # Unpack results into separate columns
    variants["FEATURE_TYPE"], variants["FUNC_SITE_DESC"], variants["IN_FUNC_SITE"] = zip(*results)
    
    # Print row count and preview
    print(f"Final annotated dataset: {len(variants):,} rows.\n")
    print("Example annotations:\n")
    print(variants[["Hugo_Symbol", "HGVSp_Short", "Protein_position", 
                     "IN_FUNC_SITE", "FEATURE_TYPE", "FUNC_SITE_DESC"]].head(10), "\n")
    
    # Return df 
    return variants


def main(): 
    args = get_args() 

    # Load variants
    print("\nLoading variants file...\n")
    variants = pd.read_csv(args.input, sep="\t", low_memory=False)
    n_input = len(variants)
    print(f"Loaded {n_input:,} variants.\n")

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
            functional_sites.extend(get_functional_sites(acc))
            time.sleep(0.25)
        functional_sites_df = pd.DataFrame(functional_sites)
        functional_sites_df.to_csv(args.funcsites, sep="\t", index=False)
        print(f"Saved functional sites to: {args.funcsites}\n")

    print(f"Loaded {len(functional_sites_df):,} functional site annotations.\n")
    print(f"Example functional sites:")
    print(f"\n{functional_sites_df.head()}\n")

    # Annotate variants
    annotated = annotate_with_functional_sites(variants.copy(), functional_sites_df)

    # Verify no variants were lost
    assert len(annotated) == n_input, f"ERROR: Lost variants! Input: {n_input}, Output: {len(annotated)}"
    print(f"All {n_input:,} variants preserved.\n")

    # Summary 
    n_total = len(annotated)
    n_true = annotated["IN_FUNC_SITE"].sum()
    print(f"{n_true:,} of {n_total:,} variants overlap functional sites ({n_true/n_total:.1%}).\n")

    # Save final annotated file 
    annotated.to_csv(args.output, sep="\t", index=False) 
    print(f"Annotated variants file saved to: {args.output}\n")


if __name__ == "__main__":
    main()